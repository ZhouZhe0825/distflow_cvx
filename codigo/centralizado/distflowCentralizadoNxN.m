function [Var, opt] = distflowCentralizadoNxN(Data, Config)

%% Inicializacion
n = size(Data.Red.Branch.T,1);

G = Data.Red.Bus.v0;
NnoG = setdiff((1:n), G)';

TSalientesG = Data.Red.Branch.T;
TEntrantesG = Data.Red.Branch.T;
TnoG = Data.Red.Branch.T;
TnoG(G,:,:) = zeros(length(G),n,Config.Etapas);
TSalientesG = TSalientesG - TnoG;
TnoG(:,G,:) = zeros(n,length(G),Config.Etapas);

NoT = 1 - Data.Red.Branch.T;

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

% Eolico
indWn = find(matOverTime(Data.Gen.DFIG.I) == 1);
NindWn = setdiff((1:n),indWn);
lenWN = length(indWn);

P_mecSigWnd = zeros(1, 1, size(Data.Gen.DFIG.P_mec,3), lenWN);
P_mecWnd = zeros(1, 1, size(Data.Gen.DFIG.P_mec,3), lenWN);
n_Wnd = zeros(1, 1, size(Data.Gen.DFIG.P_mec,3), lenWN);
TgWnd = zeros(size(Data.Gen.DFIG.Tg,1), size(Data.Gen.DFIG.Tg,2), size(Data.Gen.DFIG.Tg,3), lenWN);
for wnd = 1:lenWN
    P_mecSigWnd(1,1,:,wnd) = abs(sign(Data.Gen.DFIG.P_mec(indWn(wnd),1,:)));
    P_mecWnd(1,1,:,wnd) = Data.Gen.DFIG.P_mec(indWn(wnd),1,:);
    n_Wnd(1,1,:,wnd) = Data.Gen.DFIG.n_(indWn(wnd),1,:);
    TgWnd(:,:,:,wnd) = Data.Gen.DFIG.Tg;
end

uLowApp = zeros(size(Data.Red.Bus.uLow,1),size(Data.Red.Bus.uLow,2),size(Data.Red.Bus.uLow,3), 2);
uTopApp = uLowApp;


%% Modelo programacion matematica

tic
cvx_begin

%		 for i = 1: size(Config.Centr,1)
%			 cvx_solver_settings(Config.Centr{i,1}, Config.Centr{i,2});
%		 end

	cvx_precision high
	%% Declaracion de variables y expresiones
	%Variables generales de la red
	variable P(n, n, Config.Etapas); % Potencia activa en el arco i,j
	variable Q(n, n, Config.Etapas); % Potencia reactiva en el arco i,j
	variable l(n, n, Config.Etapas); % Corriente en el arco i,j
	variable z(n, n, Config.Etapas);
	variable y(n, n, Config.Etapas) binary;
	variable w(n,1, Config.Etapas);

	variable v(n,1, Config.Etapas); % Modulo^2 de la tension
	variable cDv(n,1, Config.Etapas); % Modulo^2 de la tension
	variable nn(n,1, Config.Etapas);
	variable nv(n,1, Config.Etapas);
	variable Tap(n,1, Config.Etapas) integer;

	variable pC(n,1, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qC(n,1, Config.Etapas); % Consumo de potencia reactiva en el nodo i

	variable pN(n,1, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qN(n,1, Config.Etapas); % Consumo de potencia reactiva en el nodo i

	variable pG(n,1, Config.Etapas); % Potencia activa generada en el nodo i
	variable qG(n,1, Config.Etapas); % Potencia reactiva generada en el nodo i

	variable qCp(n,1, Config.Etapas); % reactive power demand in i
	variable Cap(n,1, Config.Etapas) integer;

	variable pGTras(n,1, Config.Etapas);
	variable qGTras(n,1, Config.Etapas);
	variable cqGTras(n,1, Config.Etapas);

	variable pCApp(n,1, Config.Etapas, 2); % real power demand in i
	variable qCApp(n,1, Config.Etapas, 2); % real power demand in i
	variable pCn(n,1, Config.Etapas, 2); % real power demand in i
	variable pCClRes(n,1, Config.Etapas); % real power demand in i
	variable qCClRes(n,1, Config.Etapas); % real power demand in i
    
	variable pWi(n,1 ,Config.Etapas);
	variable qWi(n,1 ,Config.Etapas);
    
	expression lQoL(n, n, Config.Etapas,3);
	expression lNorm(n, n, Config.Etapas,3);
	expression vExpr(n,n,Config.Etapas);
	expression vApp(n, 1, Config.Etapas, 2);
	expression CapDif(n,1,Config.Etapas);
	expression TapDif(n,1,Config.Etapas);

	expression tfopt_expr(Config.Etapas); 
	expression fopt_expr; 


	%% Funcion objetivo
	CapDif(:,1,1) = Data.Red.Bus.CapIni;
	CapDif(:,1,(2:Config.Etapas)) = Cap(:,1,(2:Config.Etapas)) - Cap(:,1,(1:Config.Etapas-1));

	TapDif(:,1,1) = Data.Red.Bus.TapIni;
	TapDif(:,1,(2:Config.Etapas)) = Tap(:,1,(2:Config.Etapas)) - Tap(:,1,(1:Config.Etapas-1));

	tfopt_expr = ...
		sum(Data.Cost.piPTras .* pGTras,1) ...
		+ sum(Data.Cost.cdv .* cDv,1) ...
		+ sum(cqGTras,1) ...
		+ sum(Data.Red.cambioCap*Data.Red.Bus.indCap.*CapDif.^2,1) ...
		+ sum(Data.Red.cambioTap*Data.Red.Bus.indTap.*TapDif.^2,1) ...
		+ sum(Data.Util.betaT(:,1,:,1).*(pCn(:,1,:,1) - Data.Util.pzCnPref(:,1,:,1)).^2,1) ...
		;

	cqGTras >= - Data.Cost.piQmtras .* qGTras;
	cqGTras >= Data.Cost.piQMtras .* qGTras;


	%% Restricciones de Red
	% Restricciones de potencias por nodo
	pN == (permute(sum(Data.Red.Branch.T.*P - Data.Red.Branch.r.*l, 1),[2 1 3]) - sum(Data.Red.Branch.T.*P, 2));
	qN == (permute(sum(Data.Red.Branch.T.*Q - Data.Red.Branch.x.*l, 1),[2 1 3]) - sum(Data.Red.Branch.T.*Q, 2));

	pC == pCClRes;
	qC == qCClRes;

	pG == pGTras + pWi;
	qG == qGTras + qCp + qWi;

	pN - pC + pG == 0;
	qN - qC + qG == 0;

	% Restricciones de capacitores
	Cap >= Data.Red.Bus.CapLow;
	Cap <= Data.Red.Bus.CapTop;

	qCp >= NcpCapL.*v + NcpvL.*Cap - NcpCapLvL;
	qCp >= NcpCapT.*v + NcpvT.*Cap - NcpCapTvT;
	qCp <= NcpCapL.*v + NcpvT.*Cap - NcpCapLvT;
	qCp <= NcpCapT.*v + NcpvL.*Cap - NcpCapTvL;

	% Restricciones conica de corriente
	lQoL(:,:,:,1) = Data.Red.Branch.T.*(2*P);
	lQoL(:,:,:,2) = Data.Red.Branch.T.*(2*Q);
	lQoL(:,:,:,3) =  Data.Red.Branch.T.*(l - repmat(v, [1 n 1]));
	norms(lQoL,2,4) <= Data.Red.Branch.T.*(l + repmat(v, [1 n 1]));

	lNorm(:,:,:,1) = Data.Red.Branch.T.*(2*P);
	lNorm(:,:,:,2) = Data.Red.Branch.T.*(2*Q);
	lNorm(:,:,:,3) =  Data.Red.Branch.T.*(l - repmat(Data.Red.Bus.uTop.^2, [1 n 1]).*z);
	norms(lNorm,2,4) <= Data.Red.Branch.T.*(l + repmat(Data.Red.Bus.uTop.^2, [1 n 1]).*z);

	% Restriccion de la tension
	nn >= (1 + Tap.*Data.Red.Bus.Ntr).^2;
	nn <= (tnnTop + tnnLow).*(1 + Tap.*Data.Red.Bus.Ntr) - (tnnTop.*tnnLow);

	nv >= nn.*(Data.Red.Bus.uLow.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uLow.^2);
	nv >= nn.*(Data.Red.Bus.uTop.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uTop.^2);

	nv <= nn.*(Data.Red.Bus.uLow.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uLow.^2);
	nv <= nn.*(Data.Red.Bus.uTop.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uTop.^2);

	Tap >= Data.Red.Bus.TapLow;
	Tap <= Data.Red.Bus.TapTop;

	vExpr = (repmat(nv, [1 n 1]) - repmat(permute(v, [2 1 3]), [n 1 1])).*Data.Red.Branch.T ...
	- 2 * (Data.Red.Branch.r .* P + Data.Red.Branch.x .* Q) + (Data.Red.Branch.r.^2 + Data.Red.Branch.x.^2) .* l;

	0 >= (vExpr - Data.Red.Branch.T.*repmat(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2, [1 n 1]).*(1-z));
	0 <= (vExpr + Data.Red.Branch.T.*repmat(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2, [1 n 1]).*(1-z));

	cDv >= 0;
	cDv >= Data.Cost.m.*(v - (1+Data.Cost.delta));
	cDv >= - Data.Cost.m.*(v - (1-Data.Cost.delta));

	% Restricciones de arcos de la Red
	Data.Red.Branch.Tup.*z + Data.Red.Branch.Tup.*permute(z, [2 1 3]) == Data.Red.Branch.Tup.*y
	sum(z(:,NnoG,:).*Data.Red.Branch.T(:,NnoG,:),1) == 1; % los nodos de que no son barras estan todos conectados
	sum(z(:,G,:).*Data.Red.Branch.T(:,G,:),1) == 0; % no hay entradas hacia las barras
	sum(sum(z(G,:,:).*Data.Red.Branch.T(G,:,:))) >= 1; % al menos una barra conectada

	w(G,1,:) == 0;
	0 >= (repmat(w, [1 n]) - repmat(permute(w, [2 1 3]), [n 1 1]) + 1 - 100000*(1-z)).*TnoG;
	w >= 0;
	w <= n-1;

	%% Restricciones de generacion
	% Restricciones de Trasmision
	pGTras >= Data.Red.Bus.P0Low;
	pGTras <= Data.Red.Bus.P0Top;
	qGTras >= Data.Red.Bus.Q0Low;
	qGTras <= Data.Red.Bus.Q0Top;

	%% Restricciones de dominio
	v >= Data.Red.Bus.uLow.^2;
	v <= Data.Red.Bus.uTop.^2;
	l <= Data.Red.Branch.lTop.*z;
	z >= 0;
	z <= 1;
	y >= Data.Red.Branch.yLow;
	y <= Data.Red.Branch.yTop;

	NoT.*P == 0;
	NoT.*Q == 0;
	NoT.*l == 0;
	NoT.*z == 0;
	NoT.*y == 0;

	%% Restricciones de clientes residenciales

	for app = 1:2
		vApp(:,:,:,app) = v(:,:,:);
		uLowApp(:,:,:,app) = Data.Red.Bus.uLow(:,:,:);
		uTopApp(:,:,:,app) = Data.Red.Bus.uTop(:,:,:);
	end

	pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnLow.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn; 
	pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnTop.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

	pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnTop.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn;
	pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnLow.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

	qCApp == pCApp.*Data.Util.tgPhi;

	pCn >= Data.Util.pzCnLow;
	pCn <= Data.Util.pzCnTop;

	pCClRes >= sum(pCApp, 4);

	qCClRes >= sum(qCApp, 4);

	%% Restricciones de generadores Eolico
	if lenWN > 0

		% Variables de generadores eolico
		variable cqWi(n,1, Config.Etapas);

		variable PDFIG(5,5, Config.Etapas,lenWN);
		variable lDFIG(5,5, Config.Etapas,lenWN);
		variable pWigDFIG(5,1, Config.Etapas,lenWN);
		variable QDFIG(5,5, Config.Etapas,lenWN);
		variable qWigDFIG(5,1, Config.Etapas,lenWN);
		variable pCDFIG(5,1, Config.Etapas,lenWN);
		variable qCDFIG(5,1, Config.Etapas,lenWN);
		variable vDFIG(5,1, Config.Etapas,lenWN);
		variable sDFIG(5,1, Config.Etapas,lenWN);
		variable xiDFIG(5,1, Config.Etapas,lenWN);

		expression lQoLDFIG(5, 5, Config.Etapas, lenWN,3);

		expression sNormdfigF(1,1,Config.Etapas,lenWN,2);
		expression sNormdfigR(1,1,Config.Etapas,lenWN,2);

		expression PQNormdfigIE(1,1,Config.Etapas,lenWN,2);
		expression PQNormdfigIF(1,1,Config.Etapas,lenWN,2);

        tfopt_expr = tfopt_expr ...
			+ sum(Data.Cost.rhopWi .* pWi) ...
            + sum(cqWi) ...
		;

        cqWi >= - Data.Cost.rhomqWi .* qWi;
        cqWi >= Data.Cost.rhoMqWi .* qWi;

        pWi(indWn,1,:) == - permute(PDFIG(1,2,:,:) + PDFIG(1,3,:,:), [4 1 3 2]);
        qWi(indWn,1,:) == - permute(QDFIG(1,2,:,:) + QDFIG(1,3,:,:), [4 1 3 2]);

        pWi(NindWn,1,:) == 0;
        qWi(NindWn,1,:) == 0;
        
        
		% Modelo de Red interna
        permute(vDFIG(1,1,:,:), [4 1 3 2]) == v(indWn,1,:);

        PDFIG(1,2,:,:) == Data.Gen.DFIG.r(1,2,:,:) .* lDFIG(1,2,:,:) - pWigDFIG(2,1,:,:);
        PDFIG(1,3,:,:) == Data.Gen.DFIG.r(1,3,:,:) .* lDFIG(1,3,:,:) + pCDFIG(3,1,:,:);
        PDFIG(4,5,:,:) == Data.Gen.DFIG.r(4,5,:,:) .* lDFIG(4,5,:,:) - pWigDFIG(5,1,:,:);

        QDFIG(1,2,:,:) == Data.Gen.DFIG.x(1,2,:,:) .* lDFIG(1,2,:,:) - qWigDFIG(2,1,:,:);
        QDFIG(1,3,:,:) == Data.Gen.DFIG.x(1,3,:,:) .* lDFIG(1,3,:,:) - qCDFIG(3,1,:,:);
        QDFIG(4,5,:,:) == n_Wnd(1,1,:,:) .* Data.Gen.DFIG.x(4,5,:,:) .* lDFIG(4,5,:,:) - qWigDFIG(5,1,:,:);

        vDFIG(2,1,:,:) == vDFIG(1,1,:,:) - 2*(Data.Gen.DFIG.r(1,2,:,:) .* PDFIG(1,2,:,:) + Data.Gen.DFIG.x(1,2,:,:) .* QDFIG(1,2,:,:)) + (Data.Gen.DFIG.r(1,2,:,:).^2 + Data.Gen.DFIG.x(1,2,:,:).^2) .* lDFIG(1,2,:,:);
        vDFIG(3,1,:,:) == vDFIG(1,1,:,:) - 2*(Data.Gen.DFIG.r(1,3,:,:) .* PDFIG(1,3,:,:) + Data.Gen.DFIG.x(1,3,:,:) .* QDFIG(1,3,:,:)) + (Data.Gen.DFIG.r(1,3,:,:).^2 + Data.Gen.DFIG.x(1,3,:,:).^2) .* lDFIG(1,3,:,:);
        vDFIG(5,1,:,:) == vDFIG(4,1,:,:) - 2*(Data.Gen.DFIG.r(4,5,:,:) .* PDFIG(4,5,:,:) + n_Wnd(1,1,:,:).*Data.Gen.DFIG.x(4,5,:,:) .* QDFIG(4,5,:,:)) + (Data.Gen.DFIG.r(4,5,:,:).^2 + n_Wnd(1,1,:,:).^2 .* Data.Gen.DFIG.x(4,5,:,:).^2) .* lDFIG(4,5,:,:);

		% Corriente
        lQoLDFIG(:,:,:,:,1) = 2*PDFIG.*TgWnd;
        lQoLDFIG(:,:,:,:,2) = 2*QDFIG.*TgWnd;
        lQoLDFIG(:,:,:,:,3) = (lDFIG - repmat(vDFIG, [1 5 1 1])).*TgWnd;
        norms(lQoLDFIG,2,5) - (lDFIG + repmat(vDFIG, [1 5 1 1])).*TgWnd <= 0;

        Data.Gen.DFIG.lTop(1,3,:,:) >= lDFIG(1,3,:,:);
        Data.Gen.DFIG.lTop(4,5,:,:) >= lDFIG(4,5,:,:);

        (Data.Gen.DFIG.lTop(1,2,:,:) - lDFIG(1,2,:,:)).*P_mecSigWnd(1,1,:,:) >= 0;
        (Data.Gen.DFIG.sTop(3,1,:,:) - sDFIG(3,1,:,:)).*P_mecSigWnd(1,1,:,:) >= 0;
        (Data.Gen.DFIG.sTop(5,1,:,:) - sDFIG(5,1,:,:)).*P_mecSigWnd(1,1,:,:) >= 0;
        (Data.Gen.DFIG.xiTop(3,1,:,:) - xiDFIG(3,1,:,:)).*P_mecSigWnd(1,1,:,:) >= 0;
        (Data.Gen.DFIG.xiTop(5,1,:,:) - xiDFIG(5,1,:,:)).*P_mecSigWnd(1,1,:,:) >= 0;

        vDFIG(2,1,:,:) >= Data.Gen.DFIG.uLow(2,1,:,:).^2;
        vDFIG(2,1,:,:) <= Data.Gen.DFIG.uTop(2,1,:,:).^2;

        vDFIG(3,1,:,:) >= Data.Gen.DFIG.uLow(3,1,:,:).^2;
        vDFIG(3,1,:,:) <= Data.Gen.DFIG.uTop(3,1,:,:).^2;

        PQNormdfigIE(:,:,:,:,1) = PDFIG(1,2,:,:);
		PQNormdfigIE(:,:,:,:,2) = QDFIG(1,2,:,:);
        Data.Gen.DFIG.PQnorm(1,2,:,:) >= norms(PQNormdfigIE,2,5);

        PQNormdfigIF(:,:,:,:,1) = PDFIG(1,3,:,:);
		PQNormdfigIF(:,:,:,:,2) = QDFIG(1,3,:,:);
        Data.Gen.DFIG.PQnorm(1,3,:,:) >= norms(PQNormdfigIF,2,5);

        sNormdfigF(:,:,:,:,1) = pCDFIG(3,1,:,:);
		sNormdfigF(:,:,:,:,2) = qCDFIG(3,1,:,:);
        sDFIG(3,1,:,:) >= norms(sNormdfigF,2,5);
        xiDFIG(3,1,:,:) >= pCDFIG(3,1,:,:).^2 + qCDFIG(3,1,:,:).^2;

        sNormdfigF(:,:,:,:,1) = PDFIG(4,5,:,:);
		sNormdfigF(:,:,:,:,2) = QDFIG(4,5,:,:);
        sDFIG(5,1,:,:) >= norms(sNormdfigF,2,5);
        xiDFIG(5,1,:,:) >= PDFIG(4,5,:,:).^2 + QDFIG(4,5,:,:).^2;

        pCDFIG(3,1,:,:) == PDFIG(4,5,:,:) ...
            + (Data.Gen.DFIG.cv(5,1,:,:) .* sDFIG(5,1,:,:) + Data.Gen.DFIG.cr(5,1,:,:) .* xiDFIG(5,1,:,:)) ...
            + (Data.Gen.DFIG.cv(3,1,:,:) .* sDFIG(3,1,:,:) + Data.Gen.DFIG.cr(3,1,:,:) .* xiDFIG(3,1,:,:));
        vDFIG(5,1,:,:) == n_Wnd(1,1,:,:).^2*(Data.Gen.DFIG.N_er^2).*vDFIG(2,1,:,:);
        pWigDFIG(2,1,:,:) == P_mecWnd(1,1,:,:) ./ (1-n_Wnd(1,1,:,:)); %TODO es igual
        pWigDFIG(5,1,:,:) == -n_Wnd(1,1,:,:).*pWigDFIG(2,1,:,:);
        qWigDFIG(5,1,:,:) == n_Wnd(1,1,:,:).*(qWigDFIG(2,1,:,:));

		(1-TgWnd).*PDFIG == 0;
		(1-TgWnd).*lDFIG == 0;
		(1-TgWnd).*QDFIG == 0;
	
	else
        pWi == 0;
        qWi == 0;
    end

	fopt_expr = sum(tfopt_expr);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
Var.Red.Branch.P = P;
Var.Red.Branch.Q = Q;
Var.Red.Branch.l = l;
Var.Red.Branch.z = z;
Var.Red.Branch.y = y;
Var.Red.Bus.w = w;

Var.Red.Bus.v = v;
Var.Red.Bus.cDv = cDv;
Var.Red.Bus.nn = nn;
Var.Red.Bus.nv = nv;
Var.Red.Bus.Tap = Tap;

Var.Red.Bus.pC = pC;
Var.Red.Bus.qC = qC;

Var.Red.Bus.pN = pN;
Var.Red.Bus.qN = qN;

Var.Red.Bus.pG = pG;
Var.Red.Bus.qG = qG;

Var.Red.Bus.qCp = qCp;
Var.Red.Bus.Cap = Cap;

Var.Red.Bus.PTras = pGTras;
Var.Red.Bus.QTras = qGTras;

Var.ClRes.pCApp = pCApp;
Var.ClRes.qCApp = qCApp;
Var.ClRes.pC = pCClRes;
Var.ClRes.qC = qCClRes;

if (lenWN > 0)
	Var.Gen.Dfig.pWi = pWi;
	Var.Gen.Dfig.qWi = qWi;
	Var.Gen.Dfig.Branch.P = PDFIG;
	Var.Gen.Dfig.Branch.Q = QDFIG;
	Var.Gen.Dfig.Branch.l = lDFIG;
	Var.Gen.Dfig.Bus.v = vDFIG;
	Var.Gen.Dfig.Bus.pC = pCDFIG;
	Var.Gen.Dfig.Bus.qC = qCDFIG;
	Var.Gen.Dfig.Bus.pg = pWigDFIG;
	Var.Gen.Dfig.Bus.qg = qWigDFIG;
	Var.Gen.Dfig.Bus.s = sDFIG;
	Var.Gen.Dfig.Bus.xi = xiDFIG;
    Var.Gen.Dfig.Bus.n_Wnd = n_Wnd;
    Var.Gen.Dfig.Bus.P_mecWnd = P_mecWnd;
else
	Var.Gen.Dfig.pWi = zeros(n,Config.Etapas);
	Var.Gen.Dfig.qWi = zeros(n,Config.Etapas);
end

status = cvx_status;
opt = fopt_expr;
cvx_clear;
end
