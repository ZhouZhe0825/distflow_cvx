function [Var, opt] = distflowCentralizadoM(Data, Config)
% function [Var, opt] = distflow_reduced(Data)

VertI = VertIMat(Data.Red.Branch.T);
VertJ = VertJMat(Data.Red.Branch.T);
OutBr = OutBrMat(Data.Red.Branch.T);
InBr = InBrMat(Data.Red.Branch.T);

n = size(VertI,2);
m = size(VertI,1);

Data.Red.Branch.rm = NxNxT2MxT(VertI,VertJ,Data.Red.Branch.r);
Data.Red.Branch.xm = NxNxT2MxT(VertI,VertJ,Data.Red.Branch.x);
Data.Red.Branch.lTopm = NxNxT2MxT(VertI,VertJ,Data.Red.Branch.lTop);
Data.Red.Branch.yTopm = NxNxT2MxT(VertI,VertJ,Data.Red.Branch.yTop);
Data.Red.Branch.yLowm = NxNxT2MxT(VertI,VertJ,Data.Red.Branch.yLow);

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
G = find(Data.Gen.Tras.I == 1);
NnoG = setdiff((1:n),G);

tnnLow = (1 + Data.Red.Bus.TapLow.*Data.Red.Bus.Ntr);
tnnTop = (1 + Data.Red.Bus.TapTop.*Data.Red.Bus.Ntr);

%% Inicializacion

% Eolico
indWn = find(matOverTime(Data.Gen.DFIG.I) == 1);
NindWn = setdiff((1:n),indWn);
lenWN = length(indWn);

% P_mecSigWnd = zeros(lenWN, Config.Etapas);
% P_mecWnd = zeros(lenWN, Config.Etapas);
% n_Wnd = zeros(lenWN, Config.Etapas);

P_mecSigWnd = squeeze(abs(sign(Data.Gen.DFIG.P_mec(indWn,1,:))))';
P_mecWnd = squeeze(Data.Gen.DFIG.P_mec(indWn,1,:))';
n_Wnd = squeeze(Data.Gen.DFIG.n_(indWn,1,:))';



NcpCapL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow;
NcpCapT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop;

NcpvL = Data.Red.Bus.Ncp.*(Data.Red.Bus.uLow).^2;
NcpvT = Data.Red.Bus.Ncp.*(Data.Red.Bus.uTop).^2;

NcpCapLvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uLow).^2;
NcpCapTvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uTop).^2;
NcpCapLvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uTop).^2;
NcpCapTvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uLow).^2;



uLowApp = zeros(size(Data.Red.Bus.uLow,1),size(Data.Red.Bus.uLow,2), 2);
uTopApp = uLowApp;

tic
cvx_begin

% 		cvx_solver_settings('MSK_DPAR_MIO_MAX_TIME', 60);
%     for i = 1: size(Config.Centr,1)
%         cvx_solver_settings(Config.Centr{i,1}, Config.Centr{i,2});
%     end

    cvx_precision high
    %% Declaracion de variables y expresiones

	variable P(m, Config.Etapas); % Potencia activa en el arco i,j
	variable Q(m, Config.Etapas); % Potencia reactiva en el arco i,j
	variable l(m, Config.Etapas); % Corriente en el arco i,j
	variable pN(n, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qN(n, Config.Etapas); % Consumo de potencia reactiva en el nodo i
	variable v(n, Config.Etapas); % Modulo^2 de la tension
	variable cDv(n, Config.Etapas); % Modulo^2 de la tension
    variable nn(n, Config.Etapas); % module of square complex voltage in i
    variable nv(n, Config.Etapas); % module of square complex voltage in i
    variable w(n, Config.Etapas);
    variable Tap(n, Config.Etapas) integer;
    variable Cap(n, Config.Etapas) integer;

    variable pG(n, Config.Etapas); % Potencia activa generada en el nodo i
	variable qG(n, Config.Etapas); % Potencia reactiva generada en el nodo i

    variable pGTras(n, Config.Etapas);
    variable qGTras(n, Config.Etapas);
    
    variable qCp(n, Config.Etapas); % reactive power demand in i
    variable cQG(n, Config.Etapas);
    
	variable pC(n, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qC(n, Config.Etapas); % Consumo de potencia reactiva en el nodo i
    
    variable z(m, Config.Etapas);
    variable y(m, Config.Etapas) binary;
    
	variable pWi(n,Config.Etapas);
	variable qWi(n,Config.Etapas);
    
 	expression lQoL(m, Config.Etapas, 3);
 	expression lNorm(m, Config.Etapas, 3);
    expression vExpr(n, Config.Etapas);
    expression CapDif(n,Config.Etapas);
    expression TapDif(n,Config.Etapas);

    variable pCApp(n, Config.Etapas, 2); % real power demand in i
    variable qCApp(n, Config.Etapas, 2); % real power demand in i
    expression vApp(n, Config.Etapas, 2);
    variable pCn(n, Config.Etapas, 2); % real power demand in i
    
    
    variable pCClRes(n, Config.Etapas); % real power demand in i
    variable qCClRes(n, Config.Etapas); % real power demand in i
    
	expression tfopt_expr(Config.Etapas,1); 
	expression fopt_expr; 

    CapDif(:,1) = Data.Red.Bus.CapIni;
    CapDif(:,(2:Config.Etapas)) = Cap(:,(2:Config.Etapas)) - Cap(:,(1:Config.Etapas-1));

    TapDif(:,1) = Data.Red.Bus.TapIni;
    TapDif(:,(2:Config.Etapas)) = Tap(:,(2:Config.Etapas)) - Tap(:,(1:Config.Etapas-1));
        
	tfopt_expr = ...
        sum(Data.Cost.piPTrasm.*pGTras,1) ...
        + sum(Data.Cost.cdvm.*cDv,1) ...
        + sum(cQG,1) ...
        + sum(Data.Red.cambioCap*Data.Red.Bus.indCap.*(CapDif.^2),1) ...
        + sum(Data.Red.cambioTap*Data.Red.Bus.indTap.*(TapDif.^2),1) ...
        + sum(Data.Util.betaT(:,:,1).*((pCn(:,:,1) - Data.Util.pzCnPref(:,:,1)).^2),1) ...
        ;

% 	fopt_expr = sum(Data.Cost.piPTrasm.*pg) + sum(Data.Cost.cdvm.*cDv);


    cQG >= - Data.Cost.piQmtras .* qGTras;
    cQG >= Data.Cost.piQMtras .* qGTras;


	% Restriccion de balance de potencia consumida y generada
    pN == InBr*(P - Data.Red.Branch.rm.*l) - OutBr*P;
    qN == InBr*(Q - Data.Red.Branch.xm.*l) - OutBr*Q;

    pG == pGTras + pWi;
    qG == qGTras + qWi + qCp;
    
    pN - pC + pG == 0;
    qN - qC + qG == 0;
    
    pC == pCClRes;
    qC == qCClRes;
    
    % Restricciones de capacitores
    Cap >= Data.Red.Bus.CapLow;
    Cap <= Data.Red.Bus.CapTop;

    qCp >= NcpCapL.*v + NcpvL.*Cap - NcpCapLvL;
    qCp >= NcpCapT.*v + NcpvT.*Cap - NcpCapTvT;
    qCp <= NcpCapL.*v + NcpvT.*Cap - NcpCapLvT;
    qCp <= NcpCapT.*v + NcpvL.*Cap - NcpCapTvL;
    

	% Restriccion conica de corriente
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
    
    vExpr = VertI*nv - VertJ*v - 2 * (Data.Red.Branch.rm.*P + Data.Red.Branch.xm.*Q) + ((Data.Red.Branch.rm).^2 + (Data.Red.Branch.xm).^2) .* l;
    
    0 >= vExpr - (VertI*(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2)).*(1-z);
    0 <= vExpr + (VertI*(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2)).*(1-z);

    cDv >= 0;
    cDv >= (v - (1+Data.Cost.delta))*diag(Data.Cost.m);
    cDv >= (v - (1-Data.Cost.delta))*diag(-Data.Cost.m);

%     cDv >= Data.Cost.m.*(v - (1+Data.Cost.delta));
%     cDv >= - Data.Cost.m.*(v - (1-Data.Cost.delta));

    
    
    % Restricciones de arcos de la Red
    z(Tupm,:) + z(Tdownm,:) == y(Tupm,:);
    InBr(NnoG,:)*z == 1; % los nodos de que no son barras estan todos conectados
    InBr(G,:)*z == 0; % no hay entradas hacia las barras
    OutBr(G,:)*z >= 1; % al menos una barra conectada

	w(G,:) == 0;
	0 >= VertI*w - VertJ*w + 1 - 100000*(1-z); %% TODO TnoG
	w >= 0;
	w <= n-1;
	
    
    
	% Restricciones de generacion
    pGTras >= Data.Gen.Tras.pgLowm;
    pGTras <= Data.Gen.Tras.pgTopm;
    qGTras >= Data.Gen.Tras.qgLowm;
    qGTras <= Data.Gen.Tras.qgTopm;
	
	% Restricciones de caja
	v >= Data.Red.Bus.uLow.^2;
	v <= Data.Red.Bus.uTop.^2;
%     pC >= Data.Red.Bus.pCLowm;
%     qC >= Data.Red.Bus.qCLowm;
	l <= Data.Red.Branch.lTopm.*z;
    z >= 0;
    z <= 1;
    y >= Data.Red.Branch.yLowm;
    y <= Data.Red.Branch.yTopm;


    %% Cliente Residencial
    for app = 1:2
        vApp(:,:,app) = v(:,:);
        uLowApp(:,:,app) = Data.Red.Bus.uLow(:,:);
        uTopApp(:,:,app) = Data.Red.Bus.uTop(:,:);
    end

    pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnLow.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn; 
    pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnTop.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

    pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnTop.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn;
    pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnLow.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

    qCApp == pCApp.*Data.Util.tgPhi;

    pCn >= Data.Util.pzCnLow;
    pCn <= Data.Util.pzCnTop;

    pCClRes >= sum(pCApp,3);

    qCClRes >= sum(qCApp,3);

	%% Restricciones de generadores Eolico
	if lenWN > 0

		% Variables de generadores eolico
		variable cqWi(n, Config.Etapas);

		variable PdfigIE(lenWN,Config.Etapas);
		variable PdfigIF(lenWN,Config.Etapas);
		variable PdfigOR(lenWN,Config.Etapas);

		variable QdfigIE(lenWN,Config.Etapas);
		variable QdfigIF(lenWN,Config.Etapas);
		variable QdfigOR(lenWN,Config.Etapas);

		variable ldfigIE(lenWN,Config.Etapas);
		variable ldfigIF(lenWN,Config.Etapas);
		variable ldfigOR(lenWN,Config.Etapas);
		
		variable vdfigI(lenWN,Config.Etapas);
		variable vdfigE(lenWN,Config.Etapas);
		variable vdfigF(lenWN,Config.Etapas);
		variable vdfigO(lenWN,Config.Etapas);
		variable vdfigR(lenWN,Config.Etapas);
		
		variable pWigdfigE(lenWN,Config.Etapas);
		variable pWigdfigR(lenWN,Config.Etapas);
		
		variable qWigdfigE(lenWN,Config.Etapas);
		variable qWigdfigR(lenWN,Config.Etapas);
		
		variable pCdfigF(lenWN,Config.Etapas);
		variable qCdfigF(lenWN,Config.Etapas);
		
		variable sdfigF(lenWN,Config.Etapas);
		variable sdfigR(lenWN,Config.Etapas);

		variable xidfigF(lenWN,Config.Etapas);
		variable xidfigR(lenWN,Config.Etapas);

		expression lQoLdfigIE(lenWN,Config.Etapas,3);
		expression lQoLdfigIF(lenWN,Config.Etapas,3);
		expression lQoLdfigOR(lenWN,Config.Etapas,3);

		expression sNormdfigF(lenWN,Config.Etapas,2);
		expression sNormdfigR(lenWN,Config.Etapas,2);

		expression PQNormdfigIE(lenWN,Config.Etapas,2);
		expression PQNormdfigIF(lenWN,Config.Etapas,2);

		
        tfopt_expr = tfopt_expr ...
			+ sum(Data.Cost.rhopWi .* pWi) ...
            + sum(cqWi) ...
		;

        cqWi >= - Data.Cost.rhomqWi .* qWi;
        cqWi >= Data.Cost.rhoMqWi .* qWi;

        pWi(indWn,:) == - (PdfigIE + PdfigIF);
        qWi(indWn,:) == - (QdfigIE + QdfigIF);
        
        pWi(NindWn,:) == 0;
        qWi(NindWn,:) == 0;
        

		% Modelo de Red interna
        vdfigI == v(indWn,:);

        PdfigIE == Data.Gen.DFIG.rIE .* ldfigIE - pWigdfigE;
        PdfigIF == Data.Gen.DFIG.rIF .* ldfigIF + pCdfigF;
        PdfigOR == Data.Gen.DFIG.rOR .* ldfigOR - pWigdfigR;

        QdfigIE == Data.Gen.DFIG.xIE .* ldfigIE - qWigdfigE;
        QdfigIF == Data.Gen.DFIG.xIF .* ldfigIF - qCdfigF;
        QdfigOR == n_Wnd .* Data.Gen.DFIG.xOR .* ldfigOR - qWigdfigR;

        vdfigE == vdfigI - 2*(Data.Gen.DFIG.rIE .* PdfigIE + Data.Gen.DFIG.xIE .* QdfigIE) + (Data.Gen.DFIG.rIE.^2 + Data.Gen.DFIG.xIE.^2) .* ldfigIE;
        vdfigF == vdfigI - 2*(Data.Gen.DFIG.rIF .* PdfigIF + Data.Gen.DFIG.xIF .* QdfigIF) + (Data.Gen.DFIG.rIF.^2 + Data.Gen.DFIG.xIF.^2) .* ldfigIF;
        vdfigR == vdfigO - 2*(Data.Gen.DFIG.rOR .* PdfigOR + n_Wnd.*Data.Gen.DFIG.xOR .* QdfigOR) + (Data.Gen.DFIG.rOR.^2 + n_Wnd.^2 .* Data.Gen.DFIG.xOR.^2) .* ldfigOR;

		% Corriente
        lQoLdfigIE(:,:,1) = 2*PdfigIE;
        lQoLdfigIE(:,:,2) = 2*QdfigIE;
        lQoLdfigIE(:,:,3) = ldfigIE - vdfigI;
        norms(lQoLdfigIE,2,3) - (ldfigIE + vdfigI) <= 0;

        lQoLdfigIF(:,:,1) = 2*PdfigIF;
        lQoLdfigIF(:,:,2) = 2*QdfigIF;
        lQoLdfigIF(:,:,3) = ldfigIF - vdfigI;
        norms(lQoLdfigIF,2,3) - (ldfigIF + vdfigI) <= 0;

        lQoLdfigOR(:,:,1) = 2*PdfigOR;
        lQoLdfigOR(:,:,2) = 2*QdfigOR;
        lQoLdfigOR(:,:,3) = ldfigOR - vdfigO;
        norms(lQoLdfigOR,2,3) - (ldfigOR + vdfigO) <= 0;
		
        Data.Gen.DFIG.lTopIF >= ldfigIF;
        Data.Gen.DFIG.lTopOR >= ldfigOR;

        (Data.Gen.DFIG.lTopIE - ldfigIE).*P_mecSigWnd >= 0;
        (Data.Gen.DFIG.sTopF - sdfigF).*P_mecSigWnd >= 0;
        (Data.Gen.DFIG.sTopR - sdfigR).*P_mecSigWnd >= 0;

        vdfigE >= Data.Gen.DFIG.uLowE.^2;
        vdfigE <= Data.Gen.DFIG.uTopE.^2;

        vdfigF >= Data.Gen.DFIG.uLowF.^2;
        vdfigF <= Data.Gen.DFIG.uTopF.^2;

        PQNormdfigIE(:,:,1) = PdfigIE;
		PQNormdfigIE(:,:,2) = QdfigIE;
        Data.Gen.DFIG.PQnormIE >= norms(PQNormdfigIE,2,3);
        % Data.Gen.DFIG.PQnormIE >= norms([PdfigIE QdfigIE],2 ,2);

        PQNormdfigIF(:,:,1) = PdfigIF;
		PQNormdfigIF(:,:,2) = QdfigIF;
        Data.Gen.DFIG.PQnormIF >= norms(PQNormdfigIF,2,3);
        % Data.Gen.DFIG.PQnormIE >= norms([PdfigIF QdfigIF],2 ,2);


        sNormdfigF(:,:,1) = pCdfigF;
		sNormdfigF(:,:,2) = qCdfigF;
        sdfigF >= norms(sNormdfigF,2,3);
		% sdfigF >= norms([pCdfigF qCdfigF],2 ,2);
        xidfigF >= pCdfigF.^2 + qCdfigF.^2;

        sNormdfigR(:,:,1) = PdfigOR;
		sNormdfigR(:,:,2) = QdfigOR;
        sdfigR >= norms(sNormdfigR,2,3);
        % sdfigR >= norms([PdfigOR QdfigOR],2 ,2);
        xidfigR >= PdfigOR.^2 + QdfigOR.^2;

        pCdfigF == PdfigOR ...
            + (Data.Gen.DFIG.cvR .* sdfigR + Data.Gen.DFIG.crR .* xidfigR) ...
            + (Data.Gen.DFIG.cvF .* sdfigF + Data.Gen.DFIG.crF .* xidfigF);
        vdfigR == n_Wnd.^2*(Data.Gen.DFIG.N_er^2).*vdfigE;
        pWigdfigE == P_mecWnd ./ (1-n_Wnd); %TODO es igual
        pWigdfigR == -n_Wnd.*pWigdfigE;
        qWigdfigR == n_Wnd.*(qWigdfigE);
	
	else
	
        pWi == 0;
        qWi == 0;

    end

    fopt_expr = sum(tfopt_expr);
    minimize fopt_expr

    
cvx_end
cvx_status
toc

%% Reconstruccion de la solucion
% pasaje a NxNxT

	Var.Red.Branch.P = MxT2NxNxT(VertI,VertJ,P);
	Var.Red.Branch.Q = MxT2NxNxT(VertI,VertJ,Q);
	Var.Red.Branch.l = MxT2NxNxT(VertI,VertJ,l);
	Var.Red.Bus.v = permute(full(v), [1 3 2]);
	Var.Red.Bus.cDv = permute(full(cDv), [1 3 2]);
	Var.Red.Bus.nv = permute(full(nv), [1 3 2]);
	Var.Red.Bus.nn = permute(full(nn), [1 3 2]);
	Var.Red.Bus.Tap = permute(full(Tap), [1 3 2]);
	Var.Red.Bus.Cap = permute(full(Cap), [1 3 2]);

	Var.Red.Bus.pC = permute(full(pC), [1 3 2]);
	Var.Red.Bus.qC = permute(full(qC), [1 3 2]);

	Var.Red.Bus.pN = permute(full(pN), [1 3 2]);
	Var.Red.Bus.qN = permute(full(qN), [1 3 2]);

	Var.Red.Bus.pG = permute(full(pG), [1 3 2]);
	Var.Red.Bus.qG = permute(full(qG), [1 3 2]);
	Var.Red.Bus.qCp = permute(full(qCp), [1 3 2]);

	Var.Red.Branch.y = MxT2NxNxT(VertI,VertJ,y);
	Var.Red.Branch.z = MxT2NxNxT(VertI,VertJ,z);
	Var.Red.Bus.w = permute(full(w), [1 3 2]);

    Var.ClRes.pC = permute(full(pCClRes), [1 3 2]);
    Var.ClRes.qC = permute(full(qCClRes), [1 3 2]);

    Var.ClRes.pCApp = permute(full(pCApp), [1 4 2 3]);
    Var.ClRes.qCApp = permute(full(qCApp), [1 4 2 3]);

    Var.Red.Bus.PTras = zeros(n,1,Config.Etapas);
    Var.Red.Bus.PTras(G,1,:) = pGTras(G,:);
    Var.Red.Bus.QTras = zeros(n,1,Config.Etapas);
    Var.Red.Bus.QTras(G,1,:) = qGTras(G,:);

    if lenWN > 0
		Var.Gen.Dfig.pWi = permute(full(pWi), [1 3 2]);
		Var.Gen.Dfig.qWi = permute(full(qWi), [1 3 2]);
        
		Var.Gen.Dfig.Branch.P = zeros(5,5, Config.Etapas,lenWN);
        Var.Gen.Dfig.Branch.P(1,2,:,:) = PdfigIE;
        Var.Gen.Dfig.Branch.P(1,3,:,:) = PdfigIF;
        Var.Gen.Dfig.Branch.P(4,5,:,:) = PdfigOR;
        
		Var.Gen.Dfig.Branch.Q = zeros(5,5, Config.Etapas,lenWN);
        Var.Gen.Dfig.Branch.Q(1,2,:,:) = QdfigIE;
        Var.Gen.Dfig.Branch.Q(1,3,:,:) = QdfigIF;
        Var.Gen.Dfig.Branch.Q(4,5,:,:) = QdfigOR;
        
		Var.Gen.Dfig.Branch.l = zeros(5,5, Config.Etapas,lenWN);
        Var.Gen.Dfig.Branch.l(1,2,:,:) = ldfigIE;
        Var.Gen.Dfig.Branch.l(1,3,:,:) = ldfigIF;
        Var.Gen.Dfig.Branch.l(4,5,:,:) = ldfigOR;
        
		Var.Gen.Dfig.Bus.v = zeros(5,1, Config.Etapas,lenWN);
        Var.Gen.Dfig.Bus.v(1,1,:,:) = vdfigI;
        Var.Gen.Dfig.Bus.v(2,1,:,:) = vdfigE;
        Var.Gen.Dfig.Bus.v(3,1,:,:) = vdfigF;
        Var.Gen.Dfig.Bus.v(4,1,:,:) = vdfigO;
        Var.Gen.Dfig.Bus.v(5,1,:,:) = vdfigR;
        
		Var.Gen.Dfig.Bus.pC = zeros(5,1, Config.Etapas,lenWN);
        Var.Gen.Dfig.Bus.pC(3,1,:,:) = pCdfigF;

		Var.Gen.Dfig.Bus.qC = zeros(5,1, Config.Etapas,lenWN);
        Var.Gen.Dfig.Bus.qC(3,1,:,:) = qCdfigF;

		Var.Gen.Dfig.Bus.pg = zeros(5,1, Config.Etapas,lenWN);
        Var.Gen.Dfig.Bus.pg(2,1,:,:) = pWigdfigE;
        Var.Gen.Dfig.Bus.pg(5,1,:,:) = pWigdfigR;

		Var.Gen.Dfig.Bus.qg = zeros(5,1, Config.Etapas,lenWN);
        Var.Gen.Dfig.Bus.qg(2,1,:,:) = qWigdfigE;
        Var.Gen.Dfig.Bus.qg(5,1,:,:) = qWigdfigR;

		Var.Gen.Dfig.Bus.s = zeros(5,1, Config.Etapas,lenWN);
        Var.Gen.Dfig.Bus.s(3,1,:,:) = sdfigF;
        Var.Gen.Dfig.Bus.s(5,1,:,:) = sdfigR;

		Var.Gen.Dfig.Bus.xi = zeros(5,1, Config.Etapas,lenWN);
        Var.Gen.Dfig.Bus.xi(3,1,:,:) = xidfigF;
        Var.Gen.Dfig.Bus.xi(3,1,:,:) = xidfigR;
        
        Var.Gen.Dfig.Bus.n_Wnd = permute(n_Wnd, [4 3 2 1]);
        Var.Gen.Dfig.Bus.P_mecWnd = permute(P_mecWnd, [4 3 2 1]);
        

	else
		Var.Gen.Dfig.pWi = zeros(n,1,Config.Etapas);
		Var.Gen.Dfig.qWi = zeros(n,1,Config.Etapas);
    end
    
    
opt = fopt_expr;
cvx_clear
end