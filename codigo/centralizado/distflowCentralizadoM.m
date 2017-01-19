function [Var, opt] = distflowCentralizadoM(Data, Config)

leyenda = '---------------------------------------- distflowCentralizadoM ----------------------------------------'

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

% Fotovoltaico
Pv = find(sign(sum(abs(Data.Gen.Pv.I),2)) == 1);
nPv = length(Pv);
NotPv = find(sign(sum(abs(Data.Gen.Pv.I),2)) == 0);


%% Modelo programacion matematica

tic
cvx_begin

%     for i = 1: size(Config.Centr,1)
%         cvx_solver_settings(Config.Centr{i,1}, Config.Centr{i,2});
%     end

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

    variable pC(n, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qC(n, Config.Etapas); % Consumo de potencia reactiva en el nodo i

    variable pN(n, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qN(n, Config.Etapas); % Consumo de potencia reactiva en el nodo i

    variable pG(n, Config.Etapas); % Potencia activa generada en el nodo i
	variable qG(n, Config.Etapas); % Potencia reactiva generada en el nodo i

    variable qCp(n, Config.Etapas); % reactive power demand in i
    variable Cap(n, Config.Etapas) integer;
    
    variable pGTras(n, Config.Etapas);
    variable qGTras(n, Config.Etapas);
    variable cqGTras(n, Config.Etapas);
    
    variable pCApp(n, Config.Etapas, 2); % real power demand in i
    variable qCApp(n, Config.Etapas, 2); % real power demand in i
    variable pCn(n, Config.Etapas, 2); % real power demand in i
    variable pCClRes(n, Config.Etapas); % real power demand in i
    variable qCClRes(n, Config.Etapas); % real power demand in i

    % Variables de generadores solares
    variable pPv(n, Config.Etapas); % real power generation in i
    variable qPv(n, Config.Etapas); % reactive power generation in i
    
 	expression lQoL(m, Config.Etapas, 3);
 	expression lNorm(m, Config.Etapas, 3);
    expression vExpr(n, Config.Etapas);
    expression vApp(n, Config.Etapas, 2);
    expression CapDif(n,Config.Etapas);
    expression TapDif(n,Config.Etapas);

    expression tfopt_expr(Config.Etapas,1); 
	expression fopt_expr; 
    
    
    %% Funcion objetivo
    CapDif(:,1) = Data.Red.Bus.CapIni;
    CapDif(:,(2:Config.Etapas)) = Cap(:,(2:Config.Etapas)) - Cap(:,(1:Config.Etapas-1));

    TapDif(:,1) = Data.Red.Bus.TapIni;
    TapDif(:,(2:Config.Etapas)) = Tap(:,(2:Config.Etapas)) - Tap(:,(1:Config.Etapas-1));
        
	tfopt_expr = ...
        sum(Data.Cost.piPTras.*pGTras,1) ...
        + sum(Data.Cost.cdv.*cDv,1) ...
        + sum(cqGTras,1) ...
        + sum(Data.Red.cambioCap*Data.Red.Bus.indCap.*(CapDif.^2),1) ...
        + sum(Data.Red.cambioTap*Data.Red.Bus.indTap.*(TapDif.^2),1) ...
        + sum(Data.Util.betaT(:,:,1).*((pCn(:,:,1) - Data.Util.pzCnPref(:,:,1)).^2),1) ...
        ;

    cqGTras >= - Data.Cost.piQmtras .* qGTras;
    cqGTras >= Data.Cost.piQMtras .* qGTras;


    %% Restricciones de Red
    % Restricciones de potencias por nodo
    pN == InBr*(P - Data.Red.Branch.r.*l) - OutBr*P;
    qN == InBr*(Q - Data.Red.Branch.x.*l) - OutBr*Q;

    pC == pCClRes;
    qC == qCClRes;

    pG == pGTras + pPv;
    qG == qGTras + qCp + qPv;

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
	
    %% Restricciones de generacion
	% Restricciones de Trasmision
    pGTras >= Data.Gen.Tras.pgLow;
    pGTras <= Data.Gen.Tras.pgTop;
    qGTras >= Data.Gen.Tras.qgLow;
    qGTras <= Data.Gen.Tras.qgTop;
    
	%% Restricciones de caja
	v >= Data.Red.Bus.uLow.^2;
	v <= Data.Red.Bus.uTop.^2;
	l <= Data.Red.Branch.lTop.*z;
    z >= 0;
    z <= 1;
    y >= Data.Red.Branch.yLow;
    y <= Data.Red.Branch.yTop;

    %% Restricciones de clientes residenciales
    uLowApp = zeros(size(Data.Red.Bus.uLow,1),size(Data.Red.Bus.uLow,2), 2);
    uTopApp = uLowApp;

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

    %% Restricciones de generadores Solares
    if nPv > 0
    % Variables de generadores solares
        variable sPv(nPv, Config.Etapas);
        variable xiPv(nPv, Config.Etapas); % module of square complex current in i
        variable cqPv(n, Config.Etapas);
        expression SPvNorm(nPv, Config.Etapas,2);
		
        tfopt_expr = tfopt_expr ...
            + sum(Data.Cost.rhopPv .* pPv) ...
            + sum(cqPv) ...
        ;
        cqPv >= - Data.Cost.rhomqPv .* qPv;
        cqPv >= Data.Cost.rhoMqPv .* qPv;

        pPv(Pv,:) == (Data.Gen.Pv.pPvg(Pv,:) - (Data.Gen.Pv.cv(Pv,:).*sPv + Data.Gen.Pv.cr(Pv,:).*xiPv));

        SPvNorm(:,:,1) = pPv(Pv,:);
        SPvNorm(:,:,2) = qPv(Pv,:);
        sPv >= norms(SPvNorm,2,3);
        xiPv >= pPv(Pv,:).^2 + qPv(Pv,:).^2;

        sPv <= Data.Gen.Pv.sTop(Pv,:).*abs(sign(Data.Gen.Pv.pPvg(Pv,:)));
		
		pPv(NotPv,:) == 0;
		qPv(NotPv,:) == 0;

    else
        pPv == 0;
        qPv == 0;
    end

    fopt_expr = sum(tfopt_expr);
    minimize fopt_expr
    
cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT
Var.Red.Branch.P = MxT2NxNxT(VertI,VertJ,P);
Var.Red.Branch.Q = MxT2NxNxT(VertI,VertJ,Q);
Var.Red.Branch.l = MxT2NxNxT(VertI,VertJ,l);
Var.Red.Branch.z = MxT2NxNxT(VertI,VertJ,z);
Var.Red.Branch.y = MxT2NxNxT(VertI,VertJ,y);
Var.Red.Bus.w = permute(full(w), [1 3 2]);

Var.Red.Bus.v = permute(full(v), [1 3 2]);
Var.Red.Bus.cDv = permute(full(cDv), [1 3 2]);
Var.Red.Bus.nn = permute(full(nn), [1 3 2]);
Var.Red.Bus.nv = permute(full(nv), [1 3 2]);
Var.Red.Bus.Tap = permute(full(Tap), [1 3 2]);

Var.Red.Bus.pC = permute(full(pC), [1 3 2]);
Var.Red.Bus.qC = permute(full(qC), [1 3 2]);

Var.Red.Bus.pN = permute(full(pN), [1 3 2]);
Var.Red.Bus.qN = permute(full(qN), [1 3 2]);

Var.Red.Bus.pG = permute(full(pG), [1 3 2]);
Var.Red.Bus.qG = permute(full(qG), [1 3 2]);

Var.Red.Bus.qCp = permute(full(qCp), [1 3 2]);
Var.Red.Bus.Cap = permute(full(Cap), [1 3 2]);

Var.Red.Bus.PTras = zeros(n,1,Config.Etapas);
Var.Red.Bus.PTras(G,1,:) = pGTras(G,:);
Var.Red.Bus.QTras = zeros(n,1,Config.Etapas);
Var.Red.Bus.QTras(G,1,:) = qGTras(G,:);

Var.ClRes.pCApp = permute(full(pCApp), [1 4 2 3]);
Var.ClRes.qCApp = permute(full(qCApp), [1 4 2 3]);
Var.ClRes.pC = permute(full(pCClRes), [1 3 2]);
Var.ClRes.qC = permute(full(qCClRes), [1 3 2]);

if nPv > 0
	Var.Gen.Pv.pPv = pPv;
	Var.Gen.Pv.qPv = qPv;
    Var.Gen.Pv.s = pPv*0;
    Var.Gen.Pv.s(Pv,:) = sPv;
    Var.Gen.Pv.xi = pPv*0;
    Var.Gen.Pv.xi(Pv,:) = xiPv;

	Var.Gen.Pv.pPv = permute(full(Var.Gen.Pv.pPv), [1 3 2]);
	Var.Gen.Pv.qPv = permute(full(Var.Gen.Pv.qPv), [1 3 2]);
	Var.Gen.Pv.s = permute(full(Var.Gen.Pv.s), [1 3 2]);
	Var.Gen.Pv.xi = permute(full(Var.Gen.Pv.xi), [1 3 2]);
end

opt = fopt_expr;
cvx_clear
end