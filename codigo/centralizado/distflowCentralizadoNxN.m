function [Var, opt] = distflowCentralizadoNxN(Data, Config)

%% Inicializacion
n = size(Data.Red.Branch.T,1);
fixed = isfield(Data, 'Fixed');


G = find(Data.Gen.Tras.I == 1);
NnoG = setdiff((1:n), G)';

TSalientesG = Data.Red.Branch.T;
TEntrantesG = Data.Red.Branch.T;
TnoG = Data.Red.Branch.T;
TnoG(G,:,:) = zeros(length(G),n,Config.Etapas);
TSalientesG = TSalientesG - TnoG;
TnoG(:,G,:) = zeros(n,length(G),Config.Etapas);

NoT = 1 - Data.Red.Branch.T;

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


% Aire Acondicionado
AC = find(Data.St.AC.I == 1);
lenAC = length(AC);

% Baterias
nodSt = find(Data.St.Bat.I == 1);
lenSt = length(nodSt);

M = [];
if lenSt > 0
	M = zeros(Config.Etapas,Config.Etapas,n);
	for i = 1:n
		for et = 1: Config.Etapas-1
			Maux = Data.St.Bat.I(i)*[Data.St.Bat.m1(i,1,et) -Data.St.Bat.m2(i,1,et)/2; ...
				-Data.St.Bat.m2(i,1,et)/2 Data.St.Bat.m1(i,1,et)];
            M((et:et+1),(et:et+1),i) = M((et:et+1),(et:et+1),i) + Maux;
		end
	end
end

% Cargas No interrumpibles
nodCh = find(Data.ClNI.I == 1);
NnodCh = find(Data.ClNI.I == 0);
lenCh = length(nodCh);

% Eolico
indWn = find(matOverTime(Data.Gen.DFIG.I) == 1);
NindWn = setdiff((1:n),indWn);
lenWN = length(indWn);

P_mecSigWnd = abs(sign(Data.Gen.DFIG.P_mec));
P_mecWnd = Data.Gen.DFIG.P_mec;
n_Wnd = Data.Gen.DFIG.n_;

% Fotovoltaico
nodPv = find(matOverTime(Data.Gen.Pv.I) == 1);
lenPv = length(nodPv);


%% Modelo programacion matematica

tic
cvx_begin

	 for i = 1: size(Config.Centr,1)
		 cvx_solver_settings(Config.Centr{i,1}, Config.Centr{i,2});
	 end

	cvx_precision high
	%% Declaracion de variables y expresiones
	%Variables generales de la red
	variable P(n, n, Config.Etapas); % Potencia activa en el arco i,j
	variable Q(n, n, Config.Etapas); % Potencia reactiva en el arco i,j
	variable l(n, n, Config.Etapas); % Corriente en el arco i,j
	variable z(n, n, Config.Etapas);
	if fixed
		variable y(n, n, Config.Etapas);
	else
		variable y(n, n, Config.Etapas) binary;
	end
	variable w(n,1, Config.Etapas);

	variable v(n,1, Config.Etapas); % Modulo^2 de la tension
	variable cDv(n,1, Config.Etapas); % Modulo^2 de la tension
	variable nn(n,n, Config.Etapas);
	variable nv(n,n, Config.Etapas);
	if fixed
		variable Ntr(n,n, Config.Etapas);
	else
		variable Ntr(n,n, Config.Etapas) integer;
	end

	variable pC(n,1, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qC(n,1, Config.Etapas); % Consumo de potencia reactiva en el nodo i

	variable pN(n,1, Config.Etapas); % Consumo de potencia activa en el nodo i
	variable qN(n,1, Config.Etapas); % Consumo de potencia reactiva en el nodo i

	variable pG(n,1, Config.Etapas); % Potencia activa generada en el nodo i
	variable qG(n,1, Config.Etapas); % Potencia reactiva generada en el nodo i

	if fixed
		dual variable dPn;
		dual variable dQn;
	end


	variable qCp(n,1, Config.Etapas); % reactive power demand in i
	variable Ncp(n,1, Config.Etapas) integer;

	variable pGTras(n,1, Config.Etapas);
	variable qGTras(n,1, Config.Etapas);
	variable cqGTras(n,1, Config.Etapas);

	variable pCApp(n,1, Config.Etapas, 2); % real power demand in i
	variable qCApp(n,1, Config.Etapas, 2); % real power demand in i
	variable pCn(n,1, Config.Etapas, 2); % real power demand in i
	variable pCClRes(n,1, Config.Etapas); % real power demand in i
	variable qCClRes(n,1, Config.Etapas); % real power demand in i

	variable pStb(n,1, Config.Etapas);
	variable qStb(n,1, Config.Etapas);

	variable pCClNI(n,1, Config.Etapas);
	variable qCClNI(n,1, Config.Etapas);

	variable pWi(n,1 ,Config.Etapas);
	variable qWi(n,1 ,Config.Etapas);

	variable pPv(n, 1, Config.Etapas);
	variable qPv(n, 1, Config.Etapas);

	expression lQoL(n, n, Config.Etapas,3);
	expression lNorm(n, n, Config.Etapas,3);
	expression vExpr(n,n,Config.Etapas);
	expression vApp(n, 1, Config.Etapas, 2);
	expression NcpDif(n,1,Config.Etapas);
	expression NtrDif(n,n,Config.Etapas);

	expression tfopt_expr(Config.Etapas); 
	expression fopt_expr; 


	%% Funcion objetivo
	NcpDif(:,1,1) = Data.Red.Bus.NcpIni;
	NcpDif(:,1,(2:Config.Etapas)) = Ncp(:,1,(2:Config.Etapas)) - Ncp(:,1,(1:Config.Etapas-1));

	NtrDif(:,:,1) = Data.Red.Branch.NtrIni;
	NtrDif(:,:,(2:Config.Etapas)) = Ntr(:,:,(2:Config.Etapas)) - Ntr(:,:,(1:Config.Etapas-1));

	tfopt_expr = ...
		sum(Data.Cost.piPTras .* pGTras,1) ...
		+ sum(Data.Cost.cdv .* cDv,1) ...
		+ sum(cqGTras,1) ...
		+ sum(Data.Cost.cCap.*NcpDif.^2,1) ...
		+ sum(sum(Data.Cost.cTap.*NtrDif.^2,1),2) ...
		+ sum(sum(Data.Cost.cY.*y,1),2) ...
		+ sum(Data.Util.betaT(:,1,:,1).*(pCn(:,1,:,1) - Data.Util.pzCnPref(:,1,:,1)).^2,1) ...
		;

	cqGTras >= - Data.Cost.piQmtras .* qGTras;
	cqGTras >= Data.Cost.piQMtras .* qGTras;


	%% Restricciones de Red
	% Restricciones de potencias por nodo
	pN == (permute(sum(Data.Red.Branch.T.*P - Data.Red.Branch.r.*l, 1),[2 1 3]) - sum(Data.Red.Branch.T.*P, 2));
	qN == (permute(sum(Data.Red.Branch.T.*Q - Data.Red.Branch.x.*l, 1),[2 1 3]) - sum(Data.Red.Branch.T.*Q, 2));

	pC == pCClRes + pCClNI;
	qC == qCClRes + qCClNI;

	pG == pGTras + pStb + pWi + pPv;
	qG == qGTras + qCp + qStb + qWi + qPv;

	if fixed
		pN - pC + pG == 0 : dPn;
		qN - qC + qG == 0 : dQn;
	else
		pN - pC + pG == 0;
		qN - qC + qG == 0;
	end

	% Restricciones de capacitores
	Ncp >= Data.Red.Bus.NcpLow;
	Ncp <= Data.Red.Bus.NcpTop;

	qCp >= NcpCapL.*v + NcpvL.*Ncp - NcpCapLvL;
	qCp >= NcpCapT.*v + NcpvT.*Ncp - NcpCapTvT;
	qCp <= NcpCapL.*v + NcpvT.*Ncp - NcpCapLvT;
	qCp <= NcpCapT.*v + NcpvL.*Ncp - NcpCapTvL;

	% Restricciones conica de corriente
	lQoL(:,:,:,1) = Data.Red.Branch.T.*(2*P);
	lQoL(:,:,:,2) = Data.Red.Branch.T.*(2*Q);
	lQoL(:,:,:,3) = Data.Red.Branch.T.*(l - repmat(v, [1 n 1]));
	norms(lQoL,2,4) <= Data.Red.Branch.T.*(l + repmat(v, [1 n 1]));

	lNorm(:,:,:,1) = Data.Red.Branch.T.*(2*P);
	lNorm(:,:,:,2) = Data.Red.Branch.T.*(2*Q);
	lNorm(:,:,:,3) = Data.Red.Branch.T.*(l - repmat(Data.Red.Bus.uTop.^2, [1 n 1]).*z);
	norms(lNorm,2,4) <= Data.Red.Branch.T.*(l + repmat(Data.Red.Bus.uTop.^2, [1 n 1]).*z);

	% Restriccion de la tension
	nn == 1 + 2*Data.Red.Branch.Tap.*Ntr;

	nv >= nn.*repmat(Data.Red.Bus.uLow.^2, [1 n 1]) + (tnnLow.^2).*repmat(v, [1 n 1]) - (tnnLow.^2).*repmat(Data.Red.Bus.uLow.^2, [1 n 1]);
	nv >= nn.*repmat(Data.Red.Bus.uTop.^2, [1 n 1]) + (tnnTop.^2).*repmat(v, [1 n 1]) - (tnnTop.^2).*repmat(Data.Red.Bus.uTop.^2, [1 n 1]);

	nv <= nn.*repmat(Data.Red.Bus.uLow.^2, [1 n 1]) + (tnnTop.^2).*repmat(v, [1 n 1]) - (tnnTop.^2).*repmat(Data.Red.Bus.uLow.^2, [1 n 1]);
	nv <= nn.*repmat(Data.Red.Bus.uTop.^2, [1 n 1]) + (tnnLow.^2).*repmat(v, [1 n 1]) - (tnnLow.^2).*repmat(Data.Red.Bus.uTop.^2, [1 n 1]);

	Ntr >= Data.Red.Branch.NtrLow;
	Ntr <= Data.Red.Branch.NtrTop;

	vExpr = (nv - repmat(permute(v, [2 1 3]), [n 1 1])).*Data.Red.Branch.T ...
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
	pGTras >= Data.Gen.Tras.pgLow;
	pGTras <= Data.Gen.Tras.pgTop;
	qGTras >= Data.Gen.Tras.qgLow;
	qGTras <= Data.Gen.Tras.qgTop;

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
	uLowApp = zeros(size(Data.Red.Bus.uLow,1),size(Data.Red.Bus.uLow,2),size(Data.Red.Bus.uLow,3), 2);
	uTopApp = uLowApp;

	for app = 1:2
		vApp(:,:,:,app) = v(:,:,:);
		uLowApp(:,:,:,app) = Data.Red.Bus.uLow(:,:,:);
		uTopApp(:,:,:,app) = Data.Red.Bus.uTop(:,:,:);
        qCApp(:,:,:,app) == pCApp(:,:,:,app).*Data.Util.tgPhi(app);
	end

	pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnLow.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn; 
	pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnTop.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

	pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnTop.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn;
	pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnLow.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

	pCn >= Data.Util.pzCnLow;
	pCn <= Data.Util.pzCnTop;

	pCClRes >= sum(pCApp, 4);

	qCClRes >= sum(qCApp, 4);

	%% Restricciones de Aire Acondicionado
	if lenAC > 0
		% Temperatura Aire Acondicionado
		variable Tvar(n,1,Config.Etapas);
		expression TvarAnt(n,1,Config.Etapas);
		expression pCAC(n,1,Config.Etapas);

		TvarAnt(:,1,1) = Data.St.AC.tempIni;
		TvarAnt(:,1,(2:Config.Etapas)) = Tvar(:,1,(1:Config.Etapas-1));
		pCAC(:,:,:) = 0;
		pCAC(AC,1,:) = pCApp(AC,1,:,2);

		Tvar(AC,1,:) - TvarAnt(AC,1,:) + Data.St.AC.epsilon(AC,1,:).*(TvarAnt(AC,1,:) - Data.temp(AC,1,:))*Data.dt - Data.St.AC.eta(AC,1,:).*pCAC(AC,1,:)*Data.dt == 0;

		Tvar(AC,1,:) >= Data.St.AC.tempLow(AC,1,:);
		Tvar(AC,1,:) <= Data.St.AC.tempTop(AC,1,:);

		tfopt_expr = tfopt_expr + sum(Data.St.AC.beta.*(Tvar - Data.St.AC.tempPref).^2,1);

	end

	%% Restricciones de almacenamiento de bateria
	if lenSt > 0
		variable sStb(n, 1, Config.Etapas);
		variable pStgb(n, 1, Config.Etapas);
		variable xiStb(n, 1, Config.Etapas);
		variable EStb(n, 1, Config.Etapas);
		variable DlEStb(n, 1, Config.Etapas);

		expression StbNorm(n, 1, Config.Etapas,2);
		expression EStbAnt(n, 1, Config.Etapas);
		expression cStb(n,1,Config.Etapas);

		DlEStb <= 0;
		DlEStb <= EStb - Data.St.Bat.ETop.*Data.St.Bat.kapa;


		cStb = (Data.St.Bat.wOm + Data.St.Bat.m3.*(DlEStb.^2)); % falta termino de m2

			for i = 1:n
            cStb(i,1,Config.Etapas) = cStb(i,1,Config.Etapas) + (squeeze(pStgb(i,1,:))' * M(:,:,i) * squeeze(pStgb(i,1,:))).*Data.St.Bat.I(i);
			end

		tfopt_expr = tfopt_expr + sum(cStb,1) + ...
			sum(Data.St.Bat.beta(:,1,Config.Etapas).* ...
				((Data.St.Bat.ETop(:,1,Config.Etapas) - EStb(:,1,Config.Etapas).*Data.St.Bat.gama(:,1,Config.Etapas)).^2) + Data.St.Bat.wU(:,1,Config.Etapas),1)./Config.Etapas;

		EStbAnt(:,1,1) = Data.St.Bat.EIni(:,1,1);
		EStbAnt(:,1,(2:Config.Etapas)) = EStb(:,1,(1:Config.Etapas-1));

		pStb == pStgb - (Data.St.Bat.cv.*sStb + Data.St.Bat.cr.*xiStb);
		EStb == (1-Data.St.Bat.epsilon).*EStbAnt - Data.St.Bat.eta.*pStgb*Data.dt;

		StbNorm(:,:,:,1) = pStgb;
		StbNorm(:,:,:,2) = qStb;
		sStb >= norms(StbNorm,2,4);

		xiStb >= pStgb.^2 + qStb.^2;

		pStb <= Data.St.Bat.pgTop;
		sStb <= Data.St.Bat.sTop;
		xiStb <= Data.St.Bat.xiTop;
		EStb <= Data.St.Bat.ETop;

		pStb >= Data.St.Bat.pgLow;
		EStb >= Data.St.Bat.ELow;
	else
		pStb == 0;
		qStb == 0;
	end

	%% Restricciones de cargas no interrumpibles
	if lenCh > 0
		if fixed
			variable stCh(n,1, Config.Etapas);
			variable onCh(n,1, Config.Etapas);
		else
			variable stCh(n,1, Config.Etapas) binary;
			variable onCh(n,1, Config.Etapas) binary;
		end
		expression onNext(n,1,Config.Etapas)

		pCClNI == Data.ClNI.pC.*onCh;
		qCClNI == Data.ClNI.qC.*onCh;


		onNext(nodCh,1,(1:Config.Etapas-1)) = onCh(nodCh,1,(2:Config.Etapas));
		onNext(nodCh,1,Config.Etapas) = 0;
		stCh - onNext + onCh >= 0;

		sum(stCh(nodCh,1,:),3) == 1;
		sum(onCh(nodCh,1,:),3) == Data.ClNI.d(nodCh);
		stCh(nodCh,1,1) == onCh(nodCh,1,1);

		stCh(NnodCh,1,:) == 0;
		onCh(NnodCh,1,:) == 0;

		tfopt_expr = tfopt_expr + sum(Data.Util.betaE.*(pCClNI - Data.Util.pzCnPrefE).^2,1);
	else
		pCClNI == 0;
		qCClNI == 0;
	end

	%% Restricciones de generadores Eolico
	if lenWN > 0

		% Variables de generadores eolico
		variable cqWi(n,1, Config.Etapas);

		variable PdfigIE(n,1,Config.Etapas);
		variable PdfigIF(n,1,Config.Etapas);
		variable PdfigOR(n,1,Config.Etapas);

		variable QdfigIE(n,1,Config.Etapas);
		variable QdfigIF(n,1,Config.Etapas);
		variable QdfigOR(n,1,Config.Etapas);

		variable ldfigIE(n,1,Config.Etapas);
		variable ldfigIF(n,1,Config.Etapas);
		variable ldfigOR(n,1,Config.Etapas);
		
		variable vdfigI(n,1,Config.Etapas);
		variable vdfigE(n,1,Config.Etapas);
		variable vdfigF(n,1,Config.Etapas);
		variable vdfigO(n,1,Config.Etapas);
		variable vdfigR(n,1,Config.Etapas);
		
		variable pWigdfigE(n,1,Config.Etapas);
		variable pWigdfigR(n,1,Config.Etapas);
		
		variable qWigdfigE(n,1,Config.Etapas);
		variable qWigdfigR(n,1,Config.Etapas);
		
		variable pCdfigF(n,1,Config.Etapas);
		variable qCdfigF(n,1,Config.Etapas);

		variable sdfigF(n,1,Config.Etapas);
		variable sdfigR(n,1,Config.Etapas);

		variable xidfigF(n,1,Config.Etapas);
		variable xidfigR(n,1,Config.Etapas);

		expression lQoLdfigIE(n,1,Config.Etapas,3);
		expression lQoLdfigIF(n,1,Config.Etapas,3);
		expression lQoLdfigOR(n,1,Config.Etapas,3);

		expression sNormdfigF(n,1,Config.Etapas,2);
		expression sNormdfigR(n,1,Config.Etapas,2);

		expression PQNormdfigIE(n,1,Config.Etapas,2);
		expression PQNormdfigIF(n,1,Config.Etapas,2);


		tfopt_expr = tfopt_expr ...
			+ sum(Data.Cost.rhopWi .* pWi) ...
			+ sum(cqWi) ...
		;

		cqWi >= - Data.Cost.rhomqWi .* qWi;
		cqWi >= Data.Cost.rhoMqWi .* qWi;

		pWi(indWn,1,:) == - (PdfigIE(indWn,1,:) + PdfigIF(indWn,1,:));
		qWi(indWn,1,:) == - (QdfigIE(indWn,1,:) + QdfigIF(indWn,1,:));

		pWi(NindWn,1,:) == 0;
		qWi(NindWn,1,:) == 0;

		% Modelo de Red interna
		vdfigI(indWn,1,:) == v(indWn,1,:);
		vdfigI(NindWn,1,:) == 0;

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
		lQoLdfigIE(:,:,:,1) = 2*PdfigIE;
		lQoLdfigIE(:,:,:,2) = 2*QdfigIE;
		lQoLdfigIE(:,:,:,3) = ldfigIE - vdfigI;
		norms(lQoLdfigIE,2,4) - (ldfigIE + vdfigI) <= 0;

		lQoLdfigIF(:,:,:,1) = 2*PdfigIF;
		lQoLdfigIF(:,:,:,2) = 2*QdfigIF;
		lQoLdfigIF(:,:,:,3) = ldfigIF - vdfigI;
		norms(lQoLdfigIF,2,4) - (ldfigIF + vdfigI) <= 0;

		lQoLdfigOR(:,:,:,1) = 2*PdfigOR;
		lQoLdfigOR(:,:,:,2) = 2*QdfigOR;
		lQoLdfigOR(:,:,:,3) = ldfigOR - vdfigO;
		norms(lQoLdfigOR,2,4) - (ldfigOR + vdfigO) <= 0;
		
		Data.Gen.DFIG.lTopIF >= ldfigIF;
		Data.Gen.DFIG.lTopOR >= ldfigOR;

		(Data.Gen.DFIG.lTopIE).*P_mecSigWnd - ldfigIE >= 0;
		(Data.Gen.DFIG.sTopF).*P_mecSigWnd - sdfigF >= 0;
		(Data.Gen.DFIG.sTopR).*P_mecSigWnd - sdfigR >= 0;
		(Data.Gen.DFIG.xiTopF).*P_mecSigWnd - xidfigF >= 0;
		(Data.Gen.DFIG.xiTopR).*P_mecSigWnd - xidfigR >= 0;

		vdfigE >= Data.Gen.DFIG.uLowE.^2;
		vdfigE <= Data.Gen.DFIG.uTopE.^2;

		vdfigF >= Data.Gen.DFIG.uLowF.^2;
		vdfigF <= Data.Gen.DFIG.uTopF.^2;

		PQNormdfigIE(:,:,:,1) = PdfigIE;
		PQNormdfigIE(:,:,:,2) = QdfigIE;
		Data.Gen.DFIG.PQnormIE >= norms(PQNormdfigIE,2,4);

		PQNormdfigIF(:,:,:,1) = PdfigIF;
		PQNormdfigIF(:,:,:,2) = QdfigIF;
		Data.Gen.DFIG.PQnormIF >= norms(PQNormdfigIF,2,4);

		sNormdfigF(:,:,:,1) = pCdfigF;
		sNormdfigF(:,:,:,2) = qCdfigF;
		sdfigF >= norms(sNormdfigF,2,4);
		xidfigF >= pCdfigF.^2 + qCdfigF.^2;

		sNormdfigR(:,:,:,1) = PdfigOR;
		sNormdfigR(:,:,:,2) = QdfigOR;
		sdfigR >= norms(sNormdfigR,2,4);
		xidfigR >= PdfigOR.^2 + QdfigOR.^2;

		pCdfigF == PdfigOR ...
			+ (Data.Gen.DFIG.cvR .* sdfigR + Data.Gen.DFIG.crR .* xidfigR) ...
			+ (Data.Gen.DFIG.cvF .* sdfigF + Data.Gen.DFIG.crF .* xidfigF);
		vdfigR == n_Wnd.^2.*(Data.Gen.DFIG.N_er.^2).*vdfigE;
		pWigdfigE == P_mecWnd ./ (1-n_Wnd); %TODO es igual
		pWigdfigR == -n_Wnd.*pWigdfigE;
		qWigdfigR == n_Wnd.*(qWigdfigE);

	else
		pWi == 0;
		qWi == 0;
	end

	%% Restricciones de generadores Solares
	if lenPv > 0
	% Variables de generadores solares
		variable sPv(n, 1, Config.Etapas);
		variable xiPv(n, 1, Config.Etapas); % module of square complex current in i
		variable cqPv(n, 1, Config.Etapas);
		expression SPvNorm(n, 1, Config.Etapas,2);

		tfopt_expr = tfopt_expr ...
			+ sum(Data.Cost.rhopPv .* pPv) ...
			+ sum(cqPv) ...
		;
		cqPv >= - Data.Cost.rhomqPv .* qPv;
		cqPv >= Data.Cost.rhoMqPv .* qPv;

		pPv == Data.Gen.Pv.pPvg - (Data.Gen.Pv.cv.*sPv + Data.Gen.Pv.cr.*xiPv);

		SPvNorm(:,:,:,1) = pPv;
		SPvNorm(:,:,:,2) = qPv;
		sPv >= norms(SPvNorm,2,4);
		xiPv >= pPv.^2 + qPv.^2;

		sPv <= Data.Gen.Pv.sTop.*abs(sign(Data.Gen.Pv.pPvg));
 		xiPv <= Data.Gen.Pv.xiTop.*abs(sign(Data.Gen.Pv.pPvg));


	else
		pPv == 0;
		qPv == 0;
	end

	if fixed
		Ncp == round(Data.Fixed.Cap);
		Ntr == round(Data.Fixed.Tap); 
		stCh == round(Data.Fixed.stCh);
		onCh == round(Data.Fixed.onCh);
		y == round(Data.Fixed.y);
	end


	fopt_expr = sum(tfopt_expr);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
Var.Red.Branch.P = P;
Var.Red.Branch.Q = Q;
Var.Red.Branch.l = l;
Var.Red.Branch.z = round(z);
Var.Red.Branch.y = y;
Var.Red.Bus.w = w;

Var.Red.Bus.v = v;
Var.Red.Bus.cDv = cDv;
Var.Red.Branch.nn = nn;
Var.Red.Branch.nv = nv;
Var.Red.Branch.Ntr = Ntr;
Var.Red.Branch.Rtr = Ntr*0;

Var.Red.Bus.pC = pC;
Var.Red.Bus.qC = qC;

Var.Red.Bus.pN = pN;
Var.Red.Bus.qN = qN;

Var.Red.Bus.pG = pG;
Var.Red.Bus.qG = qG;

Var.Red.Bus.qCp = qCp;
Var.Red.Bus.Ncp = Ncp;

Var.Red.Bus.PTras = pGTras;
Var.Red.Bus.QTras = qGTras;

Var.ClRes.pCApp = pCApp;
Var.ClRes.qCApp = qCApp;
Var.ClRes.pC = pCClRes;
Var.ClRes.qC = qCClRes;

if fixed
	Var.Dual.dPn = dPn;
	Var.Dual.dQn = dQn;
end

% Aire Acondicionado
if lenAC > 0
	Var.ClRes.Tvar = Tvar;
end

% Baterias
if lenSt > 0
	Var.St.Bat.pStb = pStb;
	Var.St.Bat.pStgb = pStgb;
	Var.St.Bat.qStb = qStb;
	Var.St.Bat.sStb = sStb;
	Var.St.Bat.xiStb = xiStb;
	Var.St.Bat.EStb = EStb;
end

% Cargas No interrumpibles
if lenCh > 0
	Var.ClNI.pC = pCClNI;
	Var.ClNI.qC = qCClNI;
	Var.ClNI.on = onCh;
	Var.ClNI.start = stCh;
end

% Eolico
if (lenWN > 0)
	Var.Gen.Dfig.pWi = pWi;
	Var.Gen.Dfig.qWi = qWi;

	Var.Gen.Dfig.Branch.PIE = PdfigIE;
	Var.Gen.Dfig.Branch.PIF = PdfigIF;
	Var.Gen.Dfig.Branch.POR = PdfigOR;

	Var.Gen.Dfig.Branch.QIE = QdfigIE;
	Var.Gen.Dfig.Branch.QIF = QdfigIF;
	Var.Gen.Dfig.Branch.QOR = QdfigOR;

	Var.Gen.Dfig.Branch.lIE = ldfigIE;
	Var.Gen.Dfig.Branch.lIF = ldfigIF;
	Var.Gen.Dfig.Branch.lOR = ldfigOR;

	Var.Gen.Dfig.Bus.vI = vdfigI;
	Var.Gen.Dfig.Bus.vE = vdfigE;
	Var.Gen.Dfig.Bus.vF = vdfigF;
	Var.Gen.Dfig.Bus.vO = vdfigO;
	Var.Gen.Dfig.Bus.vR = vdfigR;

	Var.Gen.Dfig.Bus.pCF = pCdfigF;

	Var.Gen.Dfig.Bus.qCF = qCdfigF;

	Var.Gen.Dfig.Bus.pgE = pWigdfigE;
	Var.Gen.Dfig.Bus.pgR = pWigdfigR;

	Var.Gen.Dfig.Bus.qgE = qWigdfigE;
	Var.Gen.Dfig.Bus.qgR = qWigdfigR;

	Var.Gen.Dfig.Bus.sF = sdfigF;
	Var.Gen.Dfig.Bus.sR = sdfigR;

	Var.Gen.Dfig.Bus.xiF = xidfigF;
	Var.Gen.Dfig.Bus.xiR = xidfigR;

	Var.Gen.Dfig.Bus.n_Wnd = n_Wnd;
	Var.Gen.Dfig.Bus.P_mecWnd = P_mecWnd;
end

% Fotovoltaico
if lenPv > 0
	Var.Gen.Pv.pPv = pPv;
	Var.Gen.Pv.qPv = qPv;
	Var.Gen.Pv.s = sPv;
	Var.Gen.Pv.xi = xiPv;
end

status = cvx_status;
opt = fopt_expr;
cvx_clear;
end
