function [Var, opt] = distflowCentralizadoM(Data, Config)

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


% Aire Acondicionado
AC = find(Data.St.AC.I == 1);
nAC = length(AC);

% Baterias
St = find(sign(sum(abs(Data.St.Bat.I),2)) == 1);
nSt = length(St);
NotSt = find(sign(sum(abs(Data.St.Bat.I),2)) == 0);

M = [];
if nSt > 0
	M = zeros(2,2,nSt,Config.Etapas);
	for i = 1:nSt
		for et = 1: Config.Etapas
			M(:,:,i,et) = [Data.St.Bat.m1(St(i),et) -Data.St.Bat.m2(St(i),et)/2; ...
				-Data.St.Bat.m2(St(i),et)/2 Data.St.Bat.m1(St(i),et)];
		end
	end
end

% Cargas No interrumpibles
ClNI = find(Data.ClNI.I == 1);
nClNI = length(ClNI);

pLClNI = Data.Util.pzCnLowE(ClNI,:);
pTClNI = Data.Util.pzCnTopE(ClNI,:);
qLClNI = Data.Util.qzCnLowE(ClNI,:);
qTClNI = Data.Util.qzCnTopE(ClNI,:);
vLClNI = Data.Red.Bus.uLow(ClNI,:).^2;
vTClNI = Data.Red.Bus.uTop(ClNI,:).^2;

% Eolico
indWn = find(matOverTime(Data.Gen.DFIG.I) == 1);
NindWn = setdiff((1:n),indWn);
lenWN = length(indWn);

P_mecSigWnd = abs(sign(Data.Gen.DFIG.P_mec));
P_mecWnd = Data.Gen.DFIG.P_mec;
n_Wnd = Data.Gen.DFIG.n_;

% Fotovoltaico
Pv = find(sign(sum(abs(Data.Gen.Pv.I),2)) == 1);
nPv = length(Pv);
NotPv = find(sign(sum(abs(Data.Gen.Pv.I),2)) == 0);


%% Modelo programacion matematica

tic
cvx_begin

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
	
	variable pStb(n, Config.Etapas);
	variable qStb(n, Config.Etapas);
	
	variable pCClNI(n, Config.Etapas);
	variable qCClNI(n, Config.Etapas);
	
	variable pWi(n,Config.Etapas);
	variable qWi(n,Config.Etapas);

	variable pPv(n, Config.Etapas);
	variable qPv(n, Config.Etapas);

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
		+ sum(Data.Cost.cCap.*(CapDif.^2),1) ...
		+ sum(Data.Cost.cTap.*(TapDif.^2),1) ...
		+ sum(Data.Cost.cY.*y,1) ...
		+ sum(Data.Util.betaT(:,:,1).*((pCn(:,:,1) - Data.Util.pzCnPref(:,:,1)).^2),1) ...
		;

	cqGTras >= - Data.Cost.piQmtras .* qGTras;
	cqGTras >= Data.Cost.piQMtras .* qGTras;


	%% Restricciones de Red
	% Restricciones de potencias por nodo
	pN == InBr*(P - Data.Red.Branch.r.*l) - OutBr*P;
	qN == InBr*(Q - Data.Red.Branch.x.*l) - OutBr*Q;

	pC == pCClRes + pCClNI;
	qC == qCClRes + qCClNI;

	pG == pGTras + pStb + pWi + pPv;
	qG == qGTras + qCp + qStb + qWi + qPv;

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
	
	%% Restricciones de Aire Acondicionado
	if nAC > 0 
		% Temperatura Aire Acondicionado
		variable Tvar(n,Config.Etapas);
		expression TvarAnt(n,Config.Etapas);
		expression pCAC(n,Config.Etapas);
	
		TvarAnt(:,1) = Data.St.AC.tempIni;
		TvarAnt(:,(2:Config.Etapas)) = Tvar(:,(1:Config.Etapas-1));
		pCAC(:,:) = 0;
		pCAC(AC,:) = pCApp(AC,:,2);

		Tvar(AC,:) - TvarAnt(AC,:) + Data.St.AC.epsilon(AC,:).*(TvarAnt(AC,:) - Data.temp(AC,:))*Data.dt - Data.St.AC.eta(AC,:).*pCAC(AC,:)*Data.dt == 0;

		Tvar(AC,:) >= Data.St.AC.tempLow(AC,:);
		Tvar(AC,:) <= Data.St.AC.tempTop(AC,:);
		tfopt_expr = tfopt_expr + sum(Data.St.AC.beta.*(Tvar - Data.St.AC.tempPref).^2,1);

	end

	%% Restricciones de almacenamiento de bateria
	if nSt > 0
		variable sStb(nSt, Config.Etapas);
		variable pStgb(nSt, Config.Etapas);
		variable xiStb(nSt, Config.Etapas);
		variable EStb(n, Config.Etapas);
		variable DlEStb(n, Config.Etapas);
		
		expression StbNorm(nSt, Config.Etapas,2);
		expression EStbAnt(nSt, Config.Etapas);
		expression cStb(n, Config.Etapas);

		DlEStb <= 0;
		DlEStb <= EStb - Data.St.Bat.ETop.*Data.St.Bat.kapa;


		cStb = (Data.St.Bat.wOm + Data.St.Bat.m3.*(DlEStb.^2)).*Data.St.Bat.I; % falta termino de m2

		% Modelado de Config.Etapas
		for et = 1: Config.Etapas
			for i = 1:nSt
				j = St(i);
				if et == 1
					cStb(j,et) = cStb(j,et) + [pStgb(i,et) 0] * M(:,:,i,et) * [pStgb(i,et); 0];
				else
					cStb(j,et) = cStb(j,et) + [pStgb(i,et) pStgb(i,et-1)] * M(:,:,i,et) * [pStgb(i,et); pStgb(i,et-1)];
				end
			end
		end
		
		tfopt_expr = tfopt_expr + sum(cStb,1);
		
		tfopt_expr = tfopt_expr + ...
			sum(Data.St.Bat.I(:,Config.Etapas).*Data.St.Bat.beta(:,Config.Etapas).* ...
				((Data.St.Bat.ETop(:,Config.Etapas) - EStb(:,Config.Etapas).*Data.St.Bat.gama(:,Config.Etapas)).^2) + Data.St.Bat.wU(:,Config.Etapas),1)./Config.Etapas;

		EStbAnt(:,1) = Data.St.Bat.EIni(St,1);
		EStbAnt(:,(2:Config.Etapas)) = EStb(St,(1:Config.Etapas-1));

		pStb(St,:) == pStgb - (Data.St.Bat.cv(St,:).*sStb + Data.St.Bat.cr(St,:).*xiStb);
		EStb(St,:) == (1-Data.St.Bat.epsilon(St,:)).*EStbAnt - Data.St.Bat.eta(St,:).*pStgb*Data.dt;
		
		StbNorm(:,:,1) = pStgb;
		StbNorm(:,:,2) = qStb(St,:);
		sStb >= norms(StbNorm,2,3);

		xiStb >= pStgb.^2 + qStb(St,:).^2;

		pStb(St,:) <= Data.St.Bat.pgTop(St,:);
		sStb <= Data.St.Bat.sTop(St,:);
		xiStb <= Data.St.Bat.xiTop(St,:);
		EStb(St,:) <= Data.St.Bat.ETop(St,:);

		pStb(St,:) >= Data.St.Bat.pgLow(St,:);
		EStb(St,:) >= Data.St.Bat.ELow(St,:);

		pStb(NotSt,:) == 0;
		qStb(NotSt,:) == 0;
		EStb(NotSt,:) == 0;
	else
		pStb == 0;
		qStb == 0;
	end

	%% Restricciones de cargas no interrumpibles
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

		PdfigIE == Data.Gen.DFIG.rIE(indWn,:) .* ldfigIE - pWigdfigE;
		PdfigIF == Data.Gen.DFIG.rIF(indWn,:) .* ldfigIF + pCdfigF;
		PdfigOR == Data.Gen.DFIG.rOR(indWn,:) .* ldfigOR - pWigdfigR;

		QdfigIE == Data.Gen.DFIG.xIE(indWn,:) .* ldfigIE - qWigdfigE;
		QdfigIF == Data.Gen.DFIG.xIF(indWn,:) .* ldfigIF - qCdfigF;
		QdfigOR == n_Wnd(indWn,:) .* Data.Gen.DFIG.xOR(indWn,:) .* ldfigOR - qWigdfigR;

		vdfigE == vdfigI - 2*(Data.Gen.DFIG.rIE(indWn,:) .* PdfigIE + Data.Gen.DFIG.xIE(indWn,:) .* QdfigIE) + (Data.Gen.DFIG.rIE(indWn,:).^2 + Data.Gen.DFIG.xIE(indWn,:).^2) .* ldfigIE;
		vdfigF == vdfigI - 2*(Data.Gen.DFIG.rIF(indWn,:) .* PdfigIF + Data.Gen.DFIG.xIF(indWn,:) .* QdfigIF) + (Data.Gen.DFIG.rIF(indWn,:).^2 + Data.Gen.DFIG.xIF(indWn,:).^2) .* ldfigIF;
		vdfigR == vdfigO - 2*(Data.Gen.DFIG.rOR(indWn,:) .* PdfigOR + n_Wnd(indWn,:).*Data.Gen.DFIG.xOR(indWn,:) .* QdfigOR) + (Data.Gen.DFIG.rOR(indWn,:).^2 + n_Wnd(indWn,:).^2 .* Data.Gen.DFIG.xOR(indWn,:).^2) .* ldfigOR;

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
		
		Data.Gen.DFIG.lTopIF(indWn,:) >= ldfigIF;
		Data.Gen.DFIG.lTopOR(indWn,:) >= ldfigOR;

		(Data.Gen.DFIG.lTopIE(indWn,:)).*P_mecSigWnd(indWn,:) - ldfigIE >= 0;
		(Data.Gen.DFIG.sTopF(indWn,:)).* P_mecSigWnd(indWn,:) - sdfigF >= 0;
		(Data.Gen.DFIG.sTopR(indWn,:)) .* P_mecSigWnd(indWn,:) - sdfigR >= 0;
		(Data.Gen.DFIG.xiTopF(indWn,:)).*P_mecSigWnd(indWn,:) - xidfigF >= 0;
		(Data.Gen.DFIG.xiTopR(indWn,:)).*P_mecSigWnd(indWn,:) - xidfigR >= 0;

		vdfigE >= Data.Gen.DFIG.uLowE(indWn,:).^2;
		vdfigE <= Data.Gen.DFIG.uTopE(indWn,:).^2;

		vdfigF >= Data.Gen.DFIG.uLowF(indWn,:).^2;
		vdfigF <= Data.Gen.DFIG.uTopF(indWn,:).^2;

		PQNormdfigIE(:,:,1) = PdfigIE;
		PQNormdfigIE(:,:,2) = QdfigIE;
		Data.Gen.DFIG.PQnormIE(indWn,:) >= norms(PQNormdfigIE,2,3);

		PQNormdfigIF(:,:,1) = PdfigIF;
		PQNormdfigIF(:,:,2) = QdfigIF;
		Data.Gen.DFIG.PQnormIF(indWn,:) >= norms(PQNormdfigIF,2,3);

		sNormdfigF(:,:,1) = pCdfigF;
		sNormdfigF(:,:,2) = qCdfigF;
		sdfigF >= norms(sNormdfigF,2,3);
		xidfigF >= pCdfigF.^2 + qCdfigF.^2;

		sNormdfigR(:,:,1) = PdfigOR;
		sNormdfigR(:,:,2) = QdfigOR;
		sdfigR >= norms(sNormdfigR,2,3);
		xidfigR >= PdfigOR.^2 + QdfigOR.^2;

		pCdfigF == PdfigOR ...
			+ (Data.Gen.DFIG.cvR(indWn,:) .* sdfigR + Data.Gen.DFIG.crR(indWn,:) .* xidfigR) ...
			+ (Data.Gen.DFIG.cvF(indWn,:) .* sdfigF + Data.Gen.DFIG.crF(indWn,:) .* xidfigF);
		vdfigR == n_Wnd(indWn,:).^2.*(repmat(Data.Gen.DFIG.N_er(indWn,:).^2,[1 Config.Etapas])).*vdfigE;
		pWigdfigE == P_mecWnd(indWn,:) ./ (1-n_Wnd(indWn,:)); %TODO es igual
		pWigdfigR == -n_Wnd(indWn,:).*pWigdfigE;
		qWigdfigR == n_Wnd(indWn,:).*(qWigdfigE);
	
	else
		pWi == 0;
		qWi == 0;
	end

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
		xiPv <= Data.Gen.Pv.xiTop(Pv,:).*abs(sign(Data.Gen.Pv.pPvg(Pv,:)));

		
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

Var.Red.Bus.pC	 = 	pC	;
Var.Red.Bus.qC	 = 	qC	;

Var.Red.Bus.pN	 = 	pN	;
Var.Red.Bus.qN	 = 	qN	;

Var.Red.Bus.pG	 = 	pG	;
Var.Red.Bus.qG	 = 	qG	;

Var.Red.Bus.qCp	 = 	qCp	;
Var.Red.Bus.Cap	 = 	Cap	;

Var.Red.Bus.PTras	 = 	pGTras	;
Var.Red.Bus.QTras	 = 	qGTras	;

Var.ClRes.pCApp	 = 	pCApp	;
Var.ClRes.qCApp	 = 	qCApp	;
Var.ClRes.pC	 = 	pCClRes	;
Var.ClRes.qC	 = 	qCClRes	;

% Aire Acondicionado
if nAC > 0 
	Var.ClRes.Tvar	 = 	Tvar	;
end

% Baterias
if nSt > 0
	Var.St.Bat.pStb	 = pStb;
	Var.St.Bat.pStgb	 = pStb*0;
	Var.St.Bat.pStgb(St,:)	 = pStgb;
	Var.St.Bat.qStb	 = qStb;
	Var.St.Bat.sStb	 = pStb*0;
	Var.St.Bat.sStb(St,:)	 = sStb;
	Var.St.Bat.xiStb	 = pStb*0;
	Var.St.Bat.xiStb(St,:)	 = xiStb;
	Var.St.Bat.EStb	 = EStb;
end

% Cargas No interrumpibles
if nClNI > 0
	Var.ClNI.pC	 = pCClNI;
	Var.ClNI.qC	 = qCClNI;
	Var.ClNI.on	 = pCClNI*0;
	Var.ClNI.on(ClNI,:)	 = onClNI;
	Var.ClNI.start	 = pCClNI*0;
	Var.ClNI.start(ClNI,:)	 = stClNI;

end

% Eolico
if lenWN > 0

	Var.Gen.Dfig.pWi = pWi;
	Var.Gen.Dfig.qWi = qWi;

	Var.Gen.Dfig.Branch.PIE = pWi*0;
	Var.Gen.Dfig.Branch.PIE(indWn,:) = PdfigIE;
	Var.Gen.Dfig.Branch.PIF = pWi*0;
	Var.Gen.Dfig.Branch.PIF(indWn,:) = PdfigIF;
	Var.Gen.Dfig.Branch.POR = pWi*0;
	Var.Gen.Dfig.Branch.POR(indWn,:) = PdfigOR;

	Var.Gen.Dfig.Branch.QIE = pWi*0;
	Var.Gen.Dfig.Branch.QIE(indWn,:) = QdfigIE;
	Var.Gen.Dfig.Branch.QIF = pWi*0;
	Var.Gen.Dfig.Branch.QIF(indWn,:) = QdfigIF;
	Var.Gen.Dfig.Branch.QOR = pWi*0;
	Var.Gen.Dfig.Branch.QOR(indWn,:) = QdfigOR;

	Var.Gen.Dfig.Branch.lIE = pWi*0;
	Var.Gen.Dfig.Branch.lIE(indWn,:) = ldfigIE;
	Var.Gen.Dfig.Branch.lIF = pWi*0;
	Var.Gen.Dfig.Branch.lIF(indWn,:) = ldfigIF;
	Var.Gen.Dfig.Branch.lOR = pWi*0;
	Var.Gen.Dfig.Branch.lOR(indWn,:) = ldfigOR;

	Var.Gen.Dfig.Bus.vI = pWi*0;
	Var.Gen.Dfig.Bus.vI(indWn,:) = vdfigI;
	Var.Gen.Dfig.Bus.vE = pWi*0;
	Var.Gen.Dfig.Bus.vE(indWn,:) = vdfigE;
	Var.Gen.Dfig.Bus.vF = pWi*0;
	Var.Gen.Dfig.Bus.vF(indWn,:) = vdfigF;
	Var.Gen.Dfig.Bus.vO = pWi*0;
	Var.Gen.Dfig.Bus.vO(indWn,:) = vdfigO;
	Var.Gen.Dfig.Bus.vR = pWi*0;
	Var.Gen.Dfig.Bus.vR(indWn,:) = vdfigR;

	Var.Gen.Dfig.Bus.pCF = pWi*0;
	Var.Gen.Dfig.Bus.pCF(indWn,:) = pCdfigF;

	Var.Gen.Dfig.Bus.qCF = pWi*0;
	Var.Gen.Dfig.Bus.qCF(indWn,:) = qCdfigF;

	Var.Gen.Dfig.Bus.pgE = pWi*0;
	Var.Gen.Dfig.Bus.pgE(indWn,:) = pWigdfigE;
	Var.Gen.Dfig.Bus.pgR = pWi*0;
	Var.Gen.Dfig.Bus.pgR(indWn,:) = pWigdfigR;

	Var.Gen.Dfig.Bus.qgE = pWi*0;
	Var.Gen.Dfig.Bus.qgE(indWn,:) = qWigdfigE;
	Var.Gen.Dfig.Bus.qgR = pWi*0;
	Var.Gen.Dfig.Bus.qgR(indWn,:) = qWigdfigR;

	Var.Gen.Dfig.Bus.sF = pWi*0;
	Var.Gen.Dfig.Bus.sF(indWn,:) = sdfigF;
	Var.Gen.Dfig.Bus.sR = pWi*0;
	Var.Gen.Dfig.Bus.sR(indWn,:) = sdfigR;

	Var.Gen.Dfig.Bus.xiF = pWi*0;
	Var.Gen.Dfig.Bus.xiF(indWn,:) = xidfigF;
	Var.Gen.Dfig.Bus.xiR = pWi*0;
	Var.Gen.Dfig.Bus.xiR(indWn,:) = xidfigR;

	Var.Gen.Dfig.Bus.n_Wnd = n_Wnd;
	Var.Gen.Dfig.Bus.P_mecWnd = P_mecWnd;

end

% Fotovoltaico
if nPv > 0
	Var.Gen.Pv.pPv	 = pPv;
	Var.Gen.Pv.qPv	 = qPv;
	Var.Gen.Pv.s	 = pPv*0;
	Var.Gen.Pv.s(Pv,:)	 = sPv;
	Var.Gen.Pv.xi	 = pPv*0;
	Var.Gen.Pv.xi(Pv,:)	 = xiPv;
end

status	 = cvx_status;
opt	 = fopt_expr;
cvx_clear
end