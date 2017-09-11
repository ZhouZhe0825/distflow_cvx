function [Var, status] = df_InicialM(Data)
% function [Var, opt, status] = distflow_nod_eol_util(Data, rhopgSn, rhoP, rhopgWn, rhomqgSn, rhoMqgSn, rhomQ, rhoMQ, rhomqgWn, rhoMqgWn, delta, m, cdv, findOpt, fijarPgSn)
% Data.Red.Branch.T, Data.Red.Branch.r, Data.Red.Branch.x, Data.Gen.Pv.I son matrices cuadradas n * n ...
% Data.Red.Bus.alpha, Data.Red.Bus.cv, Data.Red.Bus.cr, Data.Gen.Pv.qgTop, uTop, Data.Red.Bus.uLow, Data.Red.Bus.pcLow, Data.Red.Bus.qcLow, Data.Red.Bus.qcCap, Data.Gen.Pv.sTop, Data.Gen.Pv.pgTop son vectores columnas de largo n
% Data.Red.Bus.v0, Data.Red.Bus.Q0Top, Data.Red.Bus.Q0Low, Data.Red.Bus.P0Top, Data.Red.Bus.P0Low son escalares
% modSym, findOpt son booleano

%% Inicializacion
VertI = VertIMat(Data.Red.Branch.T);
VertJ = VertJMat(Data.Red.Branch.T);
OutBr = VertI';
InBr = VertJ';

n = size(VertI,2);
m = size(VertI,1);

BigMtr = 1e6;
BigMw = 1e6;

G = find(Data.Gen.Tras.I == 1);
NnoG = setdiff((1:n),G);

Data.Red.Branch.Tup = triu(Data.Red.Branch.T);
[Tupm, Tdownm] = getSubTindM(Data.Red.Branch.T, triu(Data.Red.Branch.T));

TrafoTras = Data.Red.Branch.Itap;
TrafoTras(NnoG,:) = 0;
TrTras = getSubTindM(Data.Red.Branch.T, TrafoTras);

[TrTrasZ,~] = getSubTindM(Data.Red.Branch.T, TrafoTras');

[TrNTras,~] = getSubTindM(Data.Red.Branch.T, Data.Red.Branch.Itap - TrafoTras - TrafoTras');

[TrReg,~] = getSubTindM(Data.Red.Branch.T, Data.Red.Branch.Itreg);

[NoTr,~] = getSubTindM(Data.Red.Branch.T, Data.Red.Branch.T - Data.Red.Branch.Itap - Data.Red.Branch.Itreg);

tnnLow = (1 + Data.Red.Branch.NtrLow.*Data.Red.Branch.Tap);
tnnTop = (1 + Data.Red.Branch.NtrTop.*Data.Red.Branch.Tap);

NcpCap = Data.Red.Bus.Cap.*Data.Red.Bus.NcpIni;

NcpvL = Data.Red.Bus.Cap.*(Data.Red.Bus.uLow).^2;
NcpvT = Data.Red.Bus.Cap.*(Data.Red.Bus.uTop).^2;

NcpCapvL = NcpCap.*(Data.Red.Bus.uLow).^2;
NcpCapvT = NcpCap.*(Data.Red.Bus.uTop).^2;

%% Modelo programacion matematica

cvx_begin quiet

	cvx_precision high
	%% Declaracion de variables y expresiones
	%Variables generales de la red
	variable P(m, 1); % Potencia activa en el arco i,j
	variable Q(m, 1); % Potencia reactiva en el arco i,j
	variable l(m, 1); % Corriente en el arco i,j
	variable z(m, 1);
	variable y(m, 1) binary;
	variable w(n, 1);

	variable v(n, 1); % Modulo^2 de la tension
	variable cDv(n, 1); % Modulo^2 de la tension
	variable nn(m, 1);
	variable nv(m, 1);
	variable Ntr(m, 1) integer;
	variable Rtr(m,1);

	variable pN(n, 1); % Consumo de potencia activa en el nodo i
	variable qN(n, 1); % Consumo de potencia reactiva en el nodo i
	variable pC(n, 1); % Consumo de potencia activa en el nodo i
	variable qC(n, 1); % Consumo de potencia reactiva en el nodo i
	variable pG(n, 1); % Consumo de potencia activa en el nodo i
	variable qG(n, 1); % Consumo de potencia reactiva en el nodo i

	variable pGTras(n, 1);
	variable qGTras(n, 1);

	variable qCp(n, 1); % reactive power demand in i
	variable Ncp(n, 1) integer;

 	expression lQoL(m, 3);
 	expression lNorm(m, 3);
	expression vExpr(n, 1);

	expression tfopt_expr(1,1); 
	expression fopt_expr; 


	%% Funcion objetivo
	tfopt_expr = ...
		sum(Data.Cost.cdv.*cDv,1) ...
		+ sum(Data.Cost.cY.*y,1) ...
		;

	%% Restricciones de Red
	% Restricciones de potencias por nodo
	pN == InBr*(P - Data.Red.Branch.r.*l) - OutBr*P;
	qN == InBr*(Q - Data.Red.Branch.x.*l) - OutBr*Q;

	pC == Data.Red.Bus.pCLow;
	qC == Data.Red.Bus.qCLow;

	pG == pGTras + Data.St.Bat.pgIni + Data.Gen.DFIG.pgIni + Data.Gen.Pv.pgIni + Data.Gen.Basic.pgIni;
	qG == qGTras + qCp + Data.St.Bat.qgIni + Data.Gen.DFIG.qgIni + Data.Gen.Pv.qgIni + Data.Gen.Basic.qgIni;

	pN - pC + pG == 0;
	qN - qC + qG == 0;

	% Restricciones de capacitores
	Ncp == Data.Red.Bus.NcpIni;

	qCp >= NcpCap.*v + NcpvL.*Ncp - NcpCapvL;
	qCp <= NcpCap.*v + NcpvT.*Ncp - NcpCapvT;

	% Restricciones conica de corriente
	lQoL(:,1) = 2 * P;
	lQoL(:,2) = 2 * Q;
	lQoL(:,3) = l - VertI*v;
	norms(lQoL,2,2) <= l + VertI*v;

	lNorm(:,1) = 2*P;
	lNorm(:,2) = 2*Q;
	lNorm(:,3) =  l - (VertI*Data.Red.Bus.uTop.^2).*z;
	norms(lNorm,2,2) <= l + (VertI*Data.Red.Bus.uTop.^2).*z;

	nv(NoTr,:) == VertI(NoTr,:)*v

	if ~isempty(TrTras)
		Tap = Data.Red.Branch.Tap(TrTras,:);
		NT = Data.Red.Branch.NtrIni(TrTras,:);

		nn(TrTras,:) >= (1 + Tap.*Ntr(TrTras,:)).^2;
		nn(TrTras,:) <= 1 + 2*Tap.*Ntr(TrTras,:) + (Tap.^2).*NT.*Ntr(TrTras,:) + BigMtr*(1-z(TrTras,:));
		nn(TrTras,:) <= 1 + 2*Tap.*Ntr(TrTras,:) + (Tap.^2).*NT.*Ntr(TrTras,:) + BigMtr*z(TrTras,:);

		nv(TrTras,:) == nn(TrTras,:).*(VertI(TrTras,:)*Data.Red.Bus.uLow.^2);

		Ntr(TrTras,:) >= -BigMtr.*(1-z(TrTras,:));
		Ntr(TrTras,:) <= BigMtr.*z(TrTras,:);
		
		Ntr(TrTras,:) == NT;
	end

	if ~isempty(TrTrasZ)
		nv(TrTrasZ,:) == (VertI(TrTrasZ,:)*v);

		Ntr(TrTrasZ,:) == 0;
	end

	if ~isempty(TrNTras)
		VIT = VertI(TrNTras,:);
		Tap = Data.Red.Branch.Tap(TrNTras,:);
		NT = Data.Red.Branch.NtrIni(TrNTras,:);
		tnLo = tnnLow(TrNTras,:);
		tnTo = tnnTop(TrNTras,:);
		

		nn(TrNTras,:) == 1 + 2*Tap.*Ntr(TrNTras,:);

		nv(TrNTras,:) >= nn(TrNTras,:).*(VIT*(Data.Red.Bus.uLow.^2)) + (tnLo.^2).*(VIT*v) - (tnLo.^2).*(VIT*Data.Red.Bus.uLow.^2);
		nv(TrNTras,:) >= nn(TrNTras,:).*(VIT*(Data.Red.Bus.uTop.^2)) + (tnTo.^2).*(VIT*v) - (tnTo.^2).*(VIT*Data.Red.Bus.uTop.^2);
		nv(TrNTras,:) <= nn(TrNTras,:).*(VIT*(Data.Red.Bus.uLow.^2)) + (tnTo.^2).*(VIT*v) - (tnTo.^2).*(VIT*Data.Red.Bus.uLow.^2);
		nv(TrNTras,:) <= nn(TrNTras,:).*(VIT*(Data.Red.Bus.uTop.^2)) + (tnLo.^2).*(VIT*v) - (tnLo.^2).*(VIT*Data.Red.Bus.uTop.^2);

		Ntr(TrNTras,:) == NT;
	end

	if ~isempty(TrReg)
		VIT = VertI(TrReg,:);
		Tap = Data.Red.Branch.Tap(TrReg,:);
		NT = Data.Red.Branch.NtrIni(TrReg,:);
		tnLo = tnnLow(TrReg,:);
		tnTo = tnnTop(TrReg,:);

		nn(TrReg,:) == 1 + 2*Tap.*Rtr(TrReg,:);

		nv(TrReg,:) >= nn(TrReg,:).*(VIT*(Data.Red.Bus.uLow.^2)) + (tnLo.^2).*(VIT*v) - (tnLo.^2).*(VIT*Data.Red.Bus.uLow.^2);
		nv(TrReg,:) >= nn(TrReg,:).*(VIT*(Data.Red.Bus.uTop.^2)) + (tnTo.^2).*(VIT*v) - (tnTo.^2).*(VIT*Data.Red.Bus.uTop.^2);
		nv(TrReg,:) <= nn(TrReg,:).*(VIT*(Data.Red.Bus.uLow.^2)) + (tnTo.^2).*(VIT*v) - (tnTo.^2).*(VIT*Data.Red.Bus.uLow.^2);
		nv(TrReg,:) <= nn(TrReg,:).*(VIT*(Data.Red.Bus.uTop.^2)) + (tnLo.^2).*(VIT*v) - (tnLo.^2).*(VIT*Data.Red.Bus.uTop.^2);

		Rtr(TrReg,:) == NT;
	end

	vExpr = nv - VertJ*v - 2 * (Data.Red.Branch.r.*P + Data.Red.Branch.x.*Q) + ((Data.Red.Branch.r).^2 + (Data.Red.Branch.x).^2) .* l;

	0 >= vExpr - (VertI*(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2)).*(1-z);
	0 <= vExpr + (VertI*(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2)).*(1-z);

	cDv >= 0;
	cDv >= Data.Cost.m.*(v - (1+Data.Cost.delta));
	cDv >= - Data.Cost.m.*(v - (1-Data.Cost.delta));

	% Restricciones de arcos de la Red
	z(Tupm,:) + z(Tdownm,:) == y(Tupm,:);
	InBr(NnoG,:)*z == 1; % los nodos de que no son barras estan todos conectados
	InBr(G,:)*z == 0; % no hay entradas hacia las barras
	OutBr(G,:)*z >= 1; % al menos una barra conectada

	w(G,:) == 0;
	0 >= VertI*w - VertJ*w + 1 - 100000*(1-z); %% TODO TnoG
	w >= 0;
	w <= n-1;

	%% Restricciones de caja
	v >= Data.Red.Bus.uLow.^2;
	v <= Data.Red.Bus.uTop.^2;
	l <= Data.Red.Branch.lTop.*z;
	z >= 0;
	z <= 1;
	y >= Data.Red.Branch.yLow;
	y <= Data.Red.Branch.yTop;

	pGTras >= Data.Gen.Tras.pgLow;
	pGTras <= Data.Gen.Tras.pgTop;
	qGTras >= Data.Gen.Tras.qgLow;
	qGTras <= Data.Gen.Tras.qgTop;
	
	fopt_expr = sum(tfopt_expr);
	minimize fopt_expr

cvx_end

%% Construccion de la estructura de solucion
% pasaje a NxNxT
Var.Red.Branch.P	 = 	P	;
Var.Red.Branch.Q	 = 	Q	;
Var.Red.Branch.l	 = 	l	;
Var.Red.Branch.z	 = 	z	;
Var.Red.Branch.y	 = 	y	;
Var.Red.Bus.w	 = 	w	;

Var.Red.Bus.v	 = 	v	;
Var.Red.Bus.cDv	 = 	cDv	;
Var.Red.Branch.nn	 = 	nn	;
Var.Red.Branch.nv	 = 	nv	;
Var.Red.Branch.Ntr	 = 	Ntr	;
Var.Red.Branch.Rtr	 = 	Rtr	;

Var.Red.Bus.pC	 = 	pC	;
Var.Red.Bus.qC	 = 	qC	;
Var.Red.Bus.pN	 = 	pN	;
Var.Red.Bus.qN	 = 	qN	;
Var.Red.Bus.pG	 = 	pG	;
Var.Red.Bus.qG	 = 	qG	;

Var.Red.Bus.qCp	 = 	qCp	;
Var.Red.Bus.Ncp	 = 	Ncp	;

Var.Red.Bus.qG = Var.Red.Bus.qCp;

status = cvx_status;
cvx_clear
end