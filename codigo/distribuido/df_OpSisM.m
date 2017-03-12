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

G = find(Data.Gen.Tras.I == 1);
NnoG = setdiff((1:n),G);

Data.Red.Branch.Tup = triu(Data.Red.Branch.T);
[rowTup, colTup, ~] = find(Data.Red.Branch.Tup == 1);
Tind = find(Data.Red.Branch.T == 1);
Tupmind = sub2ind(size(Data.Red.Branch.T), rowTup, colTup);
Tdownmind = sub2ind(size(Data.Red.Branch.T), colTup, rowTup);
Tupm = zeros(size(Tupmind));
Tdownm = zeros(size(Tdownmind));
for i=1:length(Tupmind)
	Tupm(i) = find(Tind == Tupmind(i));
	Tdownm(i) = find(Tind == Tdownmind(i));
end

tnnLow = (1 + Data.Red.Bus.TapLow.*Data.Red.Bus.Ntr);
tnnTop = (1 + Data.Red.Bus.TapTop.*Data.Red.Bus.Ntr);

NcpCapL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow;
NcpCapT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop;

NcpvL = Data.Red.Bus.Ncp.*(Data.Red.Bus.uLow).^2;
NcpvT = Data.Red.Bus.Ncp.*(Data.Red.Bus.uTop).^2;

NcpCapLvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uLow).^2;
NcpCapTvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uTop).^2;
NcpCapLvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uTop).^2;
NcpCapTvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uLow).^2;

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
	variable nn(n, Config.Etapas);
	variable nv(n, Config.Etapas);
	variable Tap(n, Config.Etapas) integer;

	variable pN(n, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qN(n, Config.Etapas); % Consumo de potencia reactiva en el nodo i

	variable qCp(n, Config.Etapas); % reactive power demand in i
	variable Cap(n, Config.Etapas) integer;

	variable pCApp(n, Config.Etapas, 2); % real power demand in i
	variable qCApp(n, Config.Etapas, 2); % real power demand in i
	variable pCn(n, Config.Etapas, 2); % real power demand in i
	variable pCClRes(n, Config.Etapas); % real power demand in i
	variable qCClRes(n, Config.Etapas); % real power demand in i

 	expression lQoL(m, Config.Etapas, 3);
 	expression lNorm(m, Config.Etapas, 3);
	expression vExpr(n, Config.Etapas);
	expression vApp(n, Config.Etapas, 2);
	expression CapDif(n,Config.Etapas);
	expression TapDif(n,Config.Etapas);

	expression tfopt_expr(Config.Etapas,1); 
	expression tfopt_virt(Config.Etapas,1);
	expression tfopt_conv(Config.Etapas,1);
	expression fopt_expr; 


	%% Funcion objetivo
	CapDif(:,1) = Data.Red.Bus.CapIni;
	CapDif(:,(2:Config.Etapas)) = Cap(:,(2:Config.Etapas)) - Cap(:,(1:Config.Etapas-1));

	TapDif(:,1) = Data.Red.Bus.TapIni;
	TapDif(:,(2:Config.Etapas)) = Tap(:,(2:Config.Etapas)) - Tap(:,(1:Config.Etapas-1));

	tfopt_expr = ...
		sum(Data.Cost.cdv.*cDv,1) ...
		+ sum(Data.Red.cambioCap*Data.Red.Bus.indCap.*(CapDif.^2),1) ...
		+ sum(Data.Red.cambioTap*Data.Red.Bus.indTap.*(TapDif.^2),1) ...
		+ sum(Data.Red.Branch.cY.*y,1) ...
		;
	tfopt_virt = sum(- DistrInfo.muT .* pN - DistrInfo.lambdaT .* qN + DistrInfo.lambdaT .* qCp,1);
	tfopt_conv = sum(norms(pN - DistrInfo.OpSis.pN,2,2) ...
					+ norms(qN - DistrInfo.OpSis.qN,2,2) ...
					+ norms(qCp - DistrInfo.OpSis.qCp,2,2),1);

	%% Restricciones de Red
	% Restricciones de potencias por nodo
	pN == InBr*(P - Data.Red.Branch.r.*l) - OutBr*P;
	qN == InBr*(Q - Data.Red.Branch.x.*l) - OutBr*Q;

	% Restricciones de capacitores
	Cap >= Data.Red.Bus.CapLow;
	Cap <= Data.Red.Bus.CapTop;

	qCp >= NcpCapL.*v + NcpvL.*Cap - NcpCapLvL;
	qCp >= NcpCapT.*v + NcpvT.*Cap - NcpCapTvT;
	qCp <= NcpCapL.*v + NcpvT.*Cap - NcpCapLvT;
	qCp <= NcpCapT.*v + NcpvL.*Cap - NcpCapTvL;

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
	nn >= (1 + Tap.*Data.Red.Bus.Ntr).^2;
	nn <= (tnnTop + tnnLow).*(1 + Tap.*Data.Red.Bus.Ntr) - (tnnTop.*tnnLow);

	nv >= nn.*(Data.Red.Bus.uLow.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uLow.^2);
	nv >= nn.*(Data.Red.Bus.uTop.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uTop.^2);

	nv <= nn.*(Data.Red.Bus.uLow.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uLow.^2);
	nv <= nn.*(Data.Red.Bus.uTop.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uTop.^2);

	Tap >= Data.Red.Bus.TapLow;
	Tap <= Data.Red.Bus.TapTop;

	vExpr = VertI*nv - VertJ*v - 2 * (Data.Red.Branch.r.*P + Data.Red.Branch.x.*Q) + ((Data.Red.Branch.r).^2 + (Data.Red.Branch.x).^2) .* l;

	0 >= vExpr - (VertI*(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2)).*(1-z);
	0 <= vExpr + (VertI*(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2)).*(1-z);

	cDv >= 0;
	cDv >= (v - (1+Data.Cost.delta))*diag(Data.Cost.m);
	cDv >= (v - (1-Data.Cost.delta))*diag(-Data.Cost.m);

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

	fopt_expr = sum(tfopt_expr + tfopt_virt + tfopt_conv);
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
Var.Red.Bus.w	 = 	w	;

Var.Red.Bus.v	 = 	v	;
Var.Red.Bus.cDv	 = 	cDv	;
Var.Red.Bus.nn	 = 	nn	;
Var.Red.Bus.nv	 = 	nv	;
Var.Red.Bus.Tap	 = 	Tap	;

Var.Red.Bus.pN	 = 	pN	;
Var.Red.Bus.qN	 = 	qN	;

Var.Red.Bus.qCp	 = 	qCp	;
Var.Red.Bus.Cap	 = 	Cap	;

Var.Red.Bus.qG = Var.Red.Bus.qCp;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_virt);
opt(1,3) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end