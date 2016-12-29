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

NcpCapL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow;
NcpCapT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop;

NcpvL = Data.Red.Bus.Ncp.*(Data.Red.Bus.uLow).^2;
NcpvT = Data.Red.Bus.Ncp.*(Data.Red.Bus.uTop).^2;

NcpCapLvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uLow).^2;
NcpCapTvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uTop).^2;
NcpCapLvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uTop).^2;
NcpCapTvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uLow).^2;



ClNI = find(Data.ClNI.I == 1);
nClNI = length(ClNI);

pLClNI = Data.Util.pzCnLowE(ClNI,:);
pTClNI = Data.Util.pzCnTopE(ClNI,:);
qLClNI = Data.Util.qzCnLowE(ClNI,:);
qTClNI = Data.Util.qzCnTopE(ClNI,:);
vLClNI = Data.Red.Bus.uLow(ClNI,:).^2;
vTClNI = Data.Red.Bus.uTop(ClNI,:).^2;

uLowApp = zeros(size(Data.Red.Bus.uLow,1),size(Data.Red.Bus.uLow,2), 2);
uTopApp = uLowApp;

iC = Data.Red.Bus.indCons;

tic
cvx_begin

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

    variable pGTras(n, Config.Etapas);
    variable qGTras(n, Config.Etapas);

    
    variable pG(n, Config.Etapas); % Potencia activa generada en el nodo i
	variable qG(n, Config.Etapas); % Potencia reactiva generada en el nodo i
    variable qCp(n, Config.Etapas); % reactive power demand in i
    variable cQG(n, Config.Etapas);
    
	variable pC(n, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qC(n, Config.Etapas); % Consumo de potencia reactiva en el nodo i
    
    variable z(m, Config.Etapas);
    variable y(m, Config.Etapas) binary;
    
    
 	expression lQoL(m, Config.Etapas, 3);
 	expression lNorm(m, Config.Etapas, 3);
    expression vExpr(n, Config.Etapas);
    expression CapDif(n,Config.Etapas);
    expression TapDif(n,Config.Etapas);

    variable pCClNI(n, Config.Etapas);
    variable qCClNI(n, Config.Etapas);

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

    pC == pCClRes + pCClNI;
    qC == qCClRes + qCClNI;

    pG == pGTras;
    qG == qGTras + qCp;

    pN - pC + pG == 0;
    qN - qC + qG == 0;
    
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
    
    %% Clientes no interrumpibles
    if nClNI > 0
		variable onClNI(nClNI,Config.Etapas) binary;
		variable stClNI(nClNI,Config.Etapas) binary;

		expression onNext(nClNI,Config.Etapas)
        expression onClNInod(n,Config.Etapas);

        onClNInod(:,:) = 0;
        onClNInod(ClNI,:) = onClNI;


        pCClNI == Data.ClNI.pC.*onClNInod;
        qCClNI == Data.ClNI.qC.*onClNInod;

        onNext(:,(1:Config.Etapas-1)) = onClNI(:,(2:Config.Etapas));
        onNext(:,Config.Etapas) = 0;
        stClNI - onNext + onClNI >= 0;

        sum(stClNI,2) == 1;
        sum(onClNI,2) == Data.ClNI.d(ClNI);
        stClNI(:,1) == onClNI(:,1);
        
        tfopt_expr = tfopt_expr + sum(Data.Util.betaE.*((pCClNI - Data.Util.pzCnPrefE).^2),1);
    else
        pCClNI == 0;
        qCClNI == 0;
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

    if nClNI > 0
        Var.ClNI.pC = pCClNI;
        Var.ClNI.qC = qCClNI;
        Var.ClNI.on = pCClNI*0;
        Var.ClNI.on(ClNI,:) = onClNI;
        Var.ClNI.start = pCClNI*0;
        Var.ClNI.start(ClNI,:) = stClNI;
        
        Var.ClNI.pC = permute(full(Var.ClNI.pC), [1 3 2]);
        Var.ClNI.qC = permute(full(Var.ClNI.qC), [1 3 2]);
        Var.ClNI.on = permute(full(Var.ClNI.on), [1 3 2]);
        Var.ClNI.start = permute(full(Var.ClNI.start), [1 3 2]);
        
    end
    
    
opt = fopt_expr;
cvx_clear
end