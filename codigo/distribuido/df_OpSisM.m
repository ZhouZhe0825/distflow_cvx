function [Var, opt, status] = df_OpSisM(Data, Config, DistrInfo)
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

NcpCapL = Data.Red.Bus.Cap.*Data.Red.Bus.NcpLow;
NcpCapT = Data.Red.Bus.Cap.*Data.Red.Bus.NcpTop;

NcpvL = Data.Red.Bus.Cap.*(Data.Red.Bus.uLow).^2;
NcpvT = Data.Red.Bus.Cap.*(Data.Red.Bus.uTop).^2;

NcpCapLvL = Data.Red.Bus.Cap.*Data.Red.Bus.NcpLow.*(Data.Red.Bus.uLow).^2;
NcpCapTvT = Data.Red.Bus.Cap.*Data.Red.Bus.NcpTop.*(Data.Red.Bus.uTop).^2;
NcpCapLvT = Data.Red.Bus.Cap.*Data.Red.Bus.NcpLow.*(Data.Red.Bus.uTop).^2;
NcpCapTvL = Data.Red.Bus.Cap.*Data.Red.Bus.NcpTop.*(Data.Red.Bus.uLow).^2;

%% Modelo programacion matematica

tic
cvx_begin quiet

	for i = 1: size(Config.Centr,1)
		cvx_solver_settings(Config.Centr{i,1}, Config.Centr{i,2});
	end

	cvx_precision high
	%% Declaracion de variables y expresiones
	%Variables generales de la red
	variable P(m, Config.Etapas); % Potencia activa en el arco i,j
	variable Q(m, Config.Etapas); % Potencia reactiva en el arco i,j
	variable l(m, Config.Etapas); % Corriente en el arco i,j
	variable z(m, Config.Etapas);
	variable y(m, Config.Etapas) binary;
	variable w(n, Config.Etapas);

	variable v(n, Config.Etapas); % Modulo^2 de la tension
	variable cDv(n, Config.Etapas); % Modulo^2 de la tension
	variable nn(m, Config.Etapas);
	variable nv(m, Config.Etapas);
	variable Ntr(m, Config.Etapas) integer;
    variable Rtr(m,Config.Etapas);

	variable pN(n, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qN(n, Config.Etapas); % Consumo de potencia reactiva en el nodo i

	variable qCp(n, Config.Etapas); % reactive power demand in i
	variable Ncp(n, Config.Etapas) integer;

 	expression lQoL(m, Config.Etapas, 3);
 	expression lNorm(m, Config.Etapas, 3);
	expression vExpr(n, Config.Etapas);
	expression NcpDif(n,Config.Etapas);
	expression NtrDif(m,Config.Etapas);
	expression yDif(m,Config.Etapas);

	expression tfopt_expr(Config.Etapas,1); 
	expression tfopt_mu(Config.Etapas,1); 
	expression tfopt_lambda(Config.Etapas,1); 
	expression tfopt_conv(Config.Etapas,1);
	expression fopt_expr; 


	%% Funcion objetivo
	NcpDif(:,1) = Ncp(:,1) - Data.Red.Bus.NcpIni;
	NcpDif(:,(2:Config.Etapas)) = Ncp(:,(2:Config.Etapas)) - Ncp(:,(1:Config.Etapas-1));

	NtrDif(:,1) = Ntr(:,1) - Data.Red.Branch.NtrIni;
	NtrDif(:,(2:Config.Etapas)) = Ntr(:,(2:Config.Etapas)) - Ntr(:,(1:Config.Etapas-1));

	yDif(:,1) = y(:,1);
	yDif(:,(2:Config.Etapas)) = y(:,(2:Config.Etapas)) - y(:,(1:Config.Etapas-1));

	tfopt_expr = ...
		sum(Data.Cost.cdv.*cDv,1) ...
		+ sum(Data.Cost.cCap.*(NcpDif.^2),1) ...
		+ sum(Data.Cost.cTap.*(NtrDif.^2),1) ...
		+ sum(Data.Cost.cY.*yDif,1) ...
		;
	tfopt_mu = sum(- DistrInfo.muT .* pN,1);
	tfopt_lambda = sum(- DistrInfo.lambdaT .* qN + DistrInfo.lambdaT .* qCp,1);
	tfopt_conv = sum(norms(pN - DistrInfo.OpSis.pN,2,2) ...
					+ norms(qN - DistrInfo.OpSis.qN,2,2) ...
					+ norms(qCp - DistrInfo.OpSis.qCp,2,2),1);

	%% Restricciones de Red
	% Restricciones de potencias por nodo
	pN == InBr*(P - Data.Red.Branch.r.*l) - OutBr*P;
	qN == InBr*(Q - Data.Red.Branch.x.*l) - OutBr*Q;

	% Restricciones de capacitores
	Ncp >= Data.Red.Bus.NcpLow;
	Ncp <= Data.Red.Bus.NcpTop;

	qCp >= NcpCapL.*v + NcpvL.*Ncp - NcpCapLvL;
	qCp >= NcpCapT.*v + NcpvT.*Ncp - NcpCapTvT;
	qCp <= NcpCapL.*v + NcpvT.*Ncp - NcpCapLvT;
	qCp <= NcpCapT.*v + NcpvL.*Ncp - NcpCapTvL;

	% Restricciones conica de corriente
	lQoL(:,:,1) = 2 * P;
	lQoL(:,:,2) = 2 * Q;
	lQoL(:,:,3) = l - VertI*v;
	norms(lQoL,2,3) <= l + VertI*v;

	lNorm(:,:,1) = 2*P;
	lNorm(:,:,2) = 2*Q;
	lNorm(:,:,3) =  l - (VertI*Data.Red.Bus.uTop.^2).*z;
	norms(lNorm,2,3) <= l + (VertI*Data.Red.Bus.uTop.^2).*z;

	% Restriccion de la tension
	% % % % % nn >= (1 + Ntr.*Data.Red.Branch.Tap).^2;
	% % % % % nn <= (tnnTop + tnnLow).*(1 + Ntr.*Data.Red.Branch.Tap) - (tnnTop.*tnnLow);

	% % % % % nv >= nn.*(VertI*(Data.Red.Bus.uLow.^2)) + (tnnLow.^2).*(VertI*v) - (tnnLow.^2).*(VertI*Data.Red.Bus.uLow.^2);
	% % % % % nv >= nn.*(VertI*(Data.Red.Bus.uTop.^2)) + (tnnTop.^2).*(VertI*v) - (tnnTop.^2).*(VertI*Data.Red.Bus.uTop.^2);

	% % % % % nv <= nn.*(VertI*(Data.Red.Bus.uLow.^2)) + (tnnTop.^2).*(VertI*v) - (tnnTop.^2).*(VertI*Data.Red.Bus.uLow.^2);
	% % % % % nv <= nn.*(VertI*(Data.Red.Bus.uTop.^2)) + (tnnLow.^2).*(VertI*v) - (tnnLow.^2).*(VertI*Data.Red.Bus.uTop.^2);

	% % % % % Ntr >= Data.Red.Branch.NtrLow;
	% % % % % Ntr <= Data.Red.Branch.NtrTop;

    nv(NoTr,:) == VertI(NoTr,:)*v

    if ~isempty(TrTras)
        Tap = Data.Red.Branch.Tap(TrTras,:);
        NT = Data.Red.Branch.NtrTop(TrTras,:);
        NL = Data.Red.Branch.NtrLow(TrTras,:);

        nn(TrTras,:) >= (1 + Tap.*Ntr(TrTras,:)).^2;
        nn(TrTras,:) <= 1 + 2*Tap.*Ntr(TrTras,:) + (Tap.^2).*NT.*Ntr(TrTras,:) + BigMtr*(1-z(TrTras,:));
        nn(TrTras,:) <= 1 + 2*Tap.*Ntr(TrTras,:) + (Tap.^2).*NL.*Ntr(TrTras,:) + BigMtr*z(TrTras,:);

        nv(TrTras,:) == nn(TrTras,:).*(VertI(TrTras,:)*Data.Red.Bus.uLow.^2);

        Ntr(TrTras,:) >= -BigMtr.*(1-z(TrTras,:));
        Ntr(TrTras,:) <= BigMtr.*z(TrTras,:);
        
        Ntr(TrTras,:) <= NT;
        Ntr(TrTras,:) >= NL;
    end

    if ~isempty(TrTrasZ)
        nv(TrTrasZ,:) == (VertI(TrTrasZ,:)*v);

        Ntr(TrTrasZ,:) == 0;
    end

    if ~isempty(TrNTras)
        VIT = VertI(TrNTras,:);
        Tap = Data.Red.Branch.Tap(TrNTras,:);
        NT = Data.Red.Branch.NtrTop(TrNTras,:);
        NL = Data.Red.Branch.NtrLow(TrNTras,:);
        tnLo = tnnLow(TrNTras,:);
        tnTo = tnnTop(TrNTras,:);
        

        nn(TrNTras,:) == 1 + 2*Tap.*Ntr(TrNTras,:);

        nv(TrNTras,:) >= nn(TrNTras,:).*(VIT*(Data.Red.Bus.uLow.^2)) + (tnLo.^2).*(VIT*v) - (tnLo.^2).*(VIT*Data.Red.Bus.uLow.^2);
        nv(TrNTras,:) >= nn(TrNTras,:).*(VIT*(Data.Red.Bus.uTop.^2)) + (tnTo.^2).*(VIT*v) - (tnTo.^2).*(VIT*Data.Red.Bus.uTop.^2);
        nv(TrNTras,:) <= nn(TrNTras,:).*(VIT*(Data.Red.Bus.uLow.^2)) + (tnTo.^2).*(VIT*v) - (tnTo.^2).*(VIT*Data.Red.Bus.uLow.^2);
        nv(TrNTras,:) <= nn(TrNTras,:).*(VIT*(Data.Red.Bus.uTop.^2)) + (tnLo.^2).*(VIT*v) - (tnLo.^2).*(VIT*Data.Red.Bus.uTop.^2);

        Ntr(TrNTras,:) >= NL;
        Ntr(TrNTras,:) <= NT;
    end

    if ~isempty(TrReg)
        VIT = VertI(TrReg,:);
        Tap = Data.Red.Branch.Tap(TrReg,:);
        NT = Data.Red.Branch.NtrTop(TrReg,:);
        NL = Data.Red.Branch.NtrLow(TrReg,:);
        tnLo = tnnLow(TrReg,:);
        tnTo = tnnTop(TrReg,:);

        nn(TrReg,:) == 1 + 2*Tap.*Rtr(TrReg,:);

        nv(TrReg,:) >= nn(TrReg,:).*(VIT*(Data.Red.Bus.uLow.^2)) + (tnLo.^2).*(VIT*v) - (tnLo.^2).*(VIT*Data.Red.Bus.uLow.^2);
        nv(TrReg,:) >= nn(TrReg,:).*(VIT*(Data.Red.Bus.uTop.^2)) + (tnTo.^2).*(VIT*v) - (tnTo.^2).*(VIT*Data.Red.Bus.uTop.^2);
        nv(TrReg,:) <= nn(TrReg,:).*(VIT*(Data.Red.Bus.uLow.^2)) + (tnTo.^2).*(VIT*v) - (tnTo.^2).*(VIT*Data.Red.Bus.uLow.^2);
        nv(TrReg,:) <= nn(TrReg,:).*(VIT*(Data.Red.Bus.uTop.^2)) + (tnLo.^2).*(VIT*v) - (tnLo.^2).*(VIT*Data.Red.Bus.uTop.^2);

        Rtr(TrReg,:) >= NL;
        Rtr(TrReg,:) <= NT;
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

	fopt_expr = sum(tfopt_expr + tfopt_mu + tfopt_lambda + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT
Var.Red.Branch.P	 = 	P	;
Var.Red.Branch.Q	 = 	Q	;
Var.Red.Branch.l	 = 	l	;
Var.Red.Branch.z	 = 	z	;
Var.Red.Branch.y	 = 	y	;
Var.Red.Branch.yDif	 = 	yDif	;
Var.Red.Bus.w	 = 	w	;

Var.Red.Bus.v	 = 	v	;
Var.Red.Bus.cDv	 = 	cDv	;
Var.Red.Branch.nn	 = 	nn	;
Var.Red.Branch.nv	 = 	nv	;
Var.Red.Branch.Ntr	 = 	Ntr	;
Var.Red.Branch.NtrDif	 = 	NtrDif	;
Var.Red.Branch.Rtr	 = 	Rtr	;

Var.Red.Bus.pN	 = 	pN	;
Var.Red.Bus.qN	 = 	qN	;

Var.Red.Bus.qCp	 = 	qCp	;
Var.Red.Bus.Ncp	 = 	Ncp	;
Var.Red.Bus.NcpDif	 = 	NcpDif	;

Var.Red.Bus.qG = Var.Red.Bus.qCp;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_mu);
opt(1,3) = sum(tfopt_lambda);
opt(1,4) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end