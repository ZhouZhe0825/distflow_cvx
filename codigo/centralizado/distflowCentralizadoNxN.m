function [Var, opt, status] = distflowCentralizadoNxN(Data, Config, findOpt, util)
% function [Var, opt, status] = distflow_nod_eol_util(Data, rhopgSn, rhoP, rhopgWn, rhomqgSn, rhoMqgSn, rhomQ, rhoMQ, rhomqgWn, rhoMqgWn, delta, m, cdv, findOpt, fijarPgSn)
% Data.Red.Branch.T, Data.Red.Branch.r, Data.Red.Branch.x, Data.Gen.Pv.I son matrices cuadradas n * n ...
% Data.Red.Bus.alpha, Data.Red.Bus.cv, Data.Red.Bus.cr, Data.Gen.Pv.qgTop, uTop, Data.Red.Bus.uLow, Data.Red.Bus.pcLow, Data.Red.Bus.qcLow, Data.Red.Bus.qcCap, Data.Gen.Pv.sTop, Data.Gen.Pv.pgTop son vectores columnas de largo n
% Data.Red.Bus.v0, Data.Red.Bus.Q0Top, Data.Red.Bus.Q0Low, Data.Red.Bus.P0Top, Data.Red.Bus.P0Low son escalares
% modSym, findOpt son booleano


% [Data] = reshapeDataL_(Data,Etapas);
% 
n = size(Data.Red.Branch.T,1);

G = Data.Red.Bus.v0;
NnoG = setdiff((1:n), G)';

TSalientesG = Data.Red.Branch.T;
TEntrantesG = Data.Red.Branch.T;
TnoG = Data.Red.Branch.T;
TnoG(G,:,:) = zeros(length(G),n,Config.Etapas);
TSalientesG = TSalientesG - TnoG;
TnoG(:,G,:) = zeros(n,length(G),Config.Etapas);
TEntrantesG = TEntrantesG - TSalientesG - TnoG;
iC = Data.Red.Bus.indCons;




FixY = Data.Red.Branch.T .* Data.Red.Branch.yTop .* Data.Red.Branch.yLow;
FreeY = Data.Red.Branch.T .* (1 - (Data.Red.Branch.yTop .* Data.Red.Branch.yLow));
NoT = 1 - Data.Red.Branch.T;

	M = zeros(2,2,n,1,Config.Etapas);
	for i = 1:n
		for et = 1: Config.Etapas
			M(:,:,i,1,et) = Data.St.Bat.I(i,1,et)*[Data.St.Bat.m1(i,1,et) -Data.St.Bat.m2(i,1,et)/2; ...
				-Data.St.Bat.m2(i,1,et)/2 Data.St.Bat.m1(i,1,et)];
		end
	end

tnnLow = (1 + Data.Red.Bus.TapLow.*Data.Red.Bus.Ntr);
tnnTop = (1 + Data.Red.Bus.TapTop.*Data.Red.Bus.Ntr);

%% Inicializacion

% Eolico
indWn = find(matOverTime(Data.Gen.DFIG.I) == 1);
lenWN = length(indWn);
nVW = intersect(NnoG, indWn);

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

nodCh = find(Data.ClNI.I == 1);
NnodCh = find(Data.ClNI.I == 0);
lenCh = length(nodCh);


nodSt = find(Data.St.Bat.I == 1);
lenSt = length(nodSt);

nodPv = find(Data.Gen.Pv.I == 1);
lenPv = length(nodPv);

% stCh1 = sym(sym('stCh', [length(Cargas), Config.Etapas]), 'real');
% onCh1 = sym(sym('onCh', [length(Cargas), Config.Etapas]), 'real');
% chg1 = sym(zeros(length(Cargas),Config.Etapas));

isSymT = isequal(matOverTime(Data.Red.Branch.T), matOverTime(Data.Red.Branch.T)');


fixed = isfield(Data, 'Fixed');

NcpCapL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow;
NcpCapT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop;

NcpvL = Data.Red.Bus.Ncp.*(Data.Red.Bus.uLow).^2;
NcpvT = Data.Red.Bus.Ncp.*(Data.Red.Bus.uTop).^2;

NcpCapLvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uLow).^2;
NcpCapTvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uTop).^2;
NcpCapLvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uTop).^2;
NcpCapTvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uLow).^2;






if length(nVW) == lenWN && isSymT
	nVW = indWn;

	%% Modelo programacion matematica

	tic
	cvx_begin

	
	% solo para mosek
% 		cvx_solver_settings('MSK_DPAR_MIO_MAX_TIME', 3600);
%         for i = 1: size(Config.Centr,1)
%             cvx_solver_settings(Config.Centr{i,1}, Config.Centr{i,2});
%         end

        cvx_precision high
        %% Declaracion de variables y expresiones
		%Variables generales de la red
        variable P(n, n, Config.Etapas); % active power from i to j
        variable Q(n, n, Config.Etapas); % reactive power from i to j
        variable l(n, n, Config.Etapas);
        
		expression pC(n,1, Config.Etapas); % real power demand in i
		expression qC(n,1, Config.Etapas); % reactive power demand in i
        
        variable pN(n,1, Config.Etapas);
        variable qN(n,1, Config.Etapas);
        
        expression pG(n,1, Config.Etapas);
        expression qG(n,1, Config.Etapas);
		
        variable qCp(n,1, Config.Etapas); % reactive power demand in i
		variable v(n,1, Config.Etapas); % module of square complex voltage in i
		variable nn(n,1, Config.Etapas); % module of square complex voltage in i
		variable nv(n,1, Config.Etapas); % module of square complex voltage in i
        
		variable w(n,1, Config.Etapas);
        variable z(n, n, Config.Etapas);

		variable PTras(n,1, Config.Etapas);
		variable QTras(n,1, Config.Etapas);

		variable cQTras(n,1, Config.Etapas);
		variable cDv(n,1, Config.Etapas);

        if fixed
            variable Cap(n,1, Config.Etapas);
            variable Tap(n,1, Config.Etapas);
            variable y(n, n, Config.Etapas);
        else
            variable Cap(n,1, Config.Etapas) integer;
            variable Tap(n,1, Config.Etapas) integer;
            variable y(n, n, Config.Etapas) binary;
        end
        
        if lenCh > 0
            if fixed
                variable stCh(n,1, Config.Etapas);
                variable onCh(n,1, Config.Etapas);
            else
                variable stCh(n,1, Config.Etapas) binary;
                variable onCh(n,1, Config.Etapas) binary;
            end
            expression onNext(n,1,Config.Etapas)
        end
        
        if fixed
            % dual variable dPn{Config.Etapas};
            % dual variable dQn{Config.Etapas};
            dual variable dPn;
            dual variable dQn;
        end
        dual variable dTvar;

        
		% Variables de generadores solares
        if lenPv > 0
            variable pPv(n, 1, Config.Etapas); % real power generation in i
            variable qPv(n, 1, Config.Etapas); % reactive power generation in i
            variable sPv(n, 1, Config.Etapas);
            variable xiPv(n, 1, Config.Etapas); % module of square complex current in i
            variable cqPv(n, 1, Config.Etapas);
        end

		if lenWN > 0

			% Variables de generadores eolico
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
		end

		if Data.Util.Func > 0
		% Variables de utilidad
			variable pCApp(n,1, Config.Etapas, 2); % real power demand in i
			variable qCApp(n,1, Config.Etapas, 2); % real power demand in i
			expression vApp(n, 1, Config.Etapas, 2);
			variable pCn(n,1, Config.Etapas, 2); % real power demand in i
			variable qoU(n,1, Config.Etapas, 2); % reactive power demand in i
			variable qzU(n,1, Config.Etapas, 2); % reactive power demand in i

		end
		variable cqWi(n,1, Config.Etapas);

		% Expresion para p y q generados que entran en la red general
		variable pWi(n,1, Config.Etapas);
		variable qWi(n,1, Config.Etapas);

		% Temperatura Aire Acondicionado
		variable Tvar(n,1,Config.Etapas);
		expression TvarAnt(n,1,Config.Etapas);
		expression pCAC(n,1,Config.Etapas);

		% Almacenamiento en baterias
        if lenSt > 0
            variable pStb(n, 1, Config.Etapas);
            variable qStb(n, 1, Config.Etapas);
            variable sStb(n, 1, Config.Etapas);
            variable pStgb(n, 1, Config.Etapas);
            variable xiStb(n, 1, Config.Etapas);
            variable EStb(n, 1, Config.Etapas);
            variable DlEStb(n, 1, Config.Etapas);
            expression EStbAnt(n, 1, Config.Etapas);
            expression cStb(n,1);
            expression tUtilStb(n,1);
        end

        % variables capMCormT(n,Config.Etapas);
        % variables capMCormL(n,Config.Etapas);
      
        

		%% Declaracion de variables y expresiones
		expression tvExpr(n,n,Config.Etapas);
		expression CapDif(n,1,Config.Etapas);
		expression TapDif(n,1,Config.Etapas);

		% i -> e  |-> r
		%  |-> f  o

		% 1 -> 2  |-> 5
		%  |-> 3  4

		% Expresion de utilidad
		expression tUtil(n,1);

		expression tfopt_expr(Config.Etapas); 
		expression fopt_expr; 

        
		expression lQoL(n, n, Config.Etapas,3);
		expression lNorm(n, n, Config.Etapas,3);
		expression lQoLDFIG(5, 5, Config.Etapas, lenWN,3);
		

		expression tUtilAC(n,1);

		fopt_expr = 0;


		%% Funcion objetivo
% % 		switch Data.Util.Func
% % 			case 1
% % %				 tUtil = -Data.Util.betaT.*(pznU - Data.Util.pznLow).^2 + Data.Util.betaT.*(Data.Util.pznLow).^2 - Data.Util.betaT.*(poU - Data.Util.poLow).^2 + Data.Util.betaT.*(Data.Util.poLow).^2;
% % 			case 2
				tUtil = -Data.Util.betaT(:,1,:,1).*(pCn(:,1,:,1) - Data.Util.pzCnPref(:,1,:,1)).^2;
				tUtilAC = -Data.St.AC.beta.*(Tvar - Data.St.AC.tempPref).^2;
% % % 				tUtil = tUtil + tUtilAC + tUtilStb;
% % % 				tUtil = tUtil + tUtilAC;
% % 					
% % 			otherwise
% % 				tUtil = zeros(n,1);
% % 		end
        
		cF = 0;

        if lenSt > 0
            tUtilStb = (-Data.St.Bat.beta(:,1,Config.Etapas).*((Data.St.Bat.ETop(:,1,Config.Etapas) - EStb(:,1,Config.Etapas)*Data.St.Bat.gama).^2) + Data.St.Bat.wU(:,1,Config.Etapas)).*Data.St.Bat.I(:,:,Config.Etapas);
            DlEStb <= 0;
            DlEStb <= EStb - Data.St.Bat.ETop*Data.St.Bat.kapa;


            cStb = (Data.St.Bat.wOm + Data.St.Bat.m3.*(DlEStb.^2)).*Data.St.Bat.I; % falta termino de m2

            % Modelado de Config.Etapas
            for et = 1: Config.Etapas

                for i = 1:n
                    if et == 1
                        cStb(i) = cStb(i) + [pStgb(i,1,et) 0] * M(:,:,i,1,et) * [pStgb(i,1,et); 0];
                    else
                        cStb(i) = cStb(i) + [pStgb(i,1,et) pStgb(i,1, et-1)] * M(:,:,i,1,et) * [pStgb(i,1,et); pStgb(i,1, et-1)];
                    end
                end
            end
        end

        CapDif(:,1,1) = Data.Red.Bus.CapIni;
        CapDif(:,1,(2:Config.Etapas)) = Cap(:,1,(2:Config.Etapas)) - Cap(:,1,(1:Config.Etapas-1));

        TapDif(:,1,1) = Data.Red.Bus.TapIni;
        TapDif(:,1,(2:Config.Etapas)) = Tap(:,1,(2:Config.Etapas)) - Tap(:,1,(1:Config.Etapas-1));
		
% % % 			+ sum(-tUtilAC) ...
		tfopt_expr = sum(...
			+ sum(Data.Cost.piPTras .* PTras) ...
			+ sum(Data.Cost.cdv .* cDv) ...
            + sum(cQTras) ...
            + sum(Data.Red.cambioCap*CapDif.^2) ...
            + sum(Data.Red.cambioTap*TapDif.^2) ...
			+ sum(-tUtil) ...
            );

        if lenCh > 0
            tfopt_expr = tfopt_expr + sum(Data.Util.betaE.*(Data.ClNI.pC.*onCh).^2);
        end
        
        if lenWN > 0
            tfopt_expr = tfopt_expr + sum(...
                + sum(Data.Cost.rhopWi .* pWi) ...
                + sum(cqWi) ...
            );
            cqWi >= - Data.Cost.rhomqWi .* qWi;
            cqWi >= Data.Cost.rhoMqWi .* qWi;
        end
        
        if lenSt > 0
            tfopt_expr = tfopt_expr + sum(...
                + sum(cStb) ...
            ) + tUtilStb;
        end
        
        if lenPv > 0
            tfopt_expr = tfopt_expr + sum(...
                sum(Data.Cost.rhopPv .* pPv) ...
                + sum(cqPv) ...
            );
            cqPv >= - Data.Cost.rhomqPv .* qPv;
            cqPv >= Data.Cost.rhoMqPv .* qPv;
            
        end
        
		cQTras >= - Data.Cost.piQmtras .* QTras;
		cQTras >= Data.Cost.piQMtras .* QTras;


		%% Restricciones de Red

        (1-Data.Gen.DFIG.I) .* pWi == 0;
        (1-Data.Gen.DFIG.I) .* qWi == 0;

        if lenWN > 0
            pWi(nVW,1,:) == - permute(PDFIG(1,2,:,:) + PDFIG(1,3,:,:), [4 1 3 2]);
            qWi(nVW,1,:) == - permute(QDFIG(1,2,:,:) + QDFIG(1,3,:,:), [4 1 3 2]);
        end
		
		
		% Restricciones generales de branches
		pN == (permute(sum(Data.Red.Branch.T.*P - Data.Red.Branch.r.*l, 1),[2 1 3]) - sum(Data.Red.Branch.T.*P, 2));
		qN == (permute(sum(Data.Red.Branch.T.*Q - Data.Red.Branch.x.*l, 1),[2 1 3]) - sum(Data.Red.Branch.T.*Q, 2));


		Cap >= Data.Red.Bus.CapLow;
		Cap <= Data.Red.Bus.CapTop;

		% Restricciones de potencias de la red
        qCp >= NcpCapL.*v + NcpvL.*Cap - NcpCapLvL;
        qCp >= NcpCapT.*v + NcpvT.*Cap - NcpCapTvT;
        qCp <= NcpCapL.*v + NcpvT.*Cap - NcpCapLvT;
        qCp <= NcpCapT.*v + NcpvL.*Cap - NcpCapTvL;


        pC = sum(pCApp, 4);
        qC = sum(qCApp, 4);
        if lenCh > 0
            pC = pC + Data.ClNI.pC.*onCh;
            qC = qC + Data.ClNI.qC.*onCh;
        end
        
        pG = PTras;
        qG = QTras + qCp;
        
        if lenWN > 0
            pG = pG + pWi;
            qG = qG + qWi;
        end
        
        if lenSt > 0
            pG = pG + pStb;
            qG = qG + qStb;
        end
        
        if lenPv > 0
            pG = pG + pPv;
            qG = qG + qPv;
        end
            
        if fixed
% 			qCp == Data.Red.Bus.Ncp.*Data.Fixed.Cap.*v;

            pN - pC + pG == 0 : dPn;
			qN - qC + qG == 0 : dQn;
		else

			pN - pC + pG == 0;
			qN - qC + qG == 0;
		end

		% Restricciones de la corriente
		lQoL(:,:,:,1) = 2*P;
		lQoL(:,:,:,2) = 2*Q;
		lQoL(:,:,:,3) =  (l - repmat(v, [1 n 1]));
		norms(lQoL,2,4) <= (l + repmat(v, [1 n 1]));

		lNorm(:,:,:,1) = 2*P;
		lNorm(:,:,:,2) = 2*Q;
		lNorm(:,:,:,3) =  (l- z.*repmat(Data.Red.Bus.uTop.^2, [1 n 1]));
		norms(lNorm,2,4) <= (l + z.*repmat(Data.Red.Bus.uTop.^2, [1 n 1]));

		l <= Data.Red.Branch.lTop.*z;

		% Restricciones de voltaje

% 		nn >= (1 + Tap.*Data.Red.Bus.Ntr).^2;
% 		nn <= (tnnTop + tnnLow).*(1 + Tap.*Data.Red.Bus.Ntr) - (tnnTop.*tnnLow);
% 
% 		nv >= nn.*(Data.Red.Bus.uLow.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uLow.^2);
% 		nv >= nn.*(Data.Red.Bus.uTop.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uTop.^2);
% 
% 		nv <= nn.*(Data.Red.Bus.uLow.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uLow.^2);
% 		nv <= nn.*(Data.Red.Bus.uTop.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uTop.^2);
% 
% 		Tap >= Data.Red.Bus.TapLow;
% 		Tap <= Data.Red.Bus.TapTop;
% 
% 		tvExpr = (repmat(nv, [1 n 1]) - repmat(permute(v, [2 1 3]), [n 1 1])).*Data.Red.Branch.T ...
% 			- 2 * (Data.Red.Branch.r .* P + Data.Red.Branch.x .* Q) + (Data.Red.Branch.r.^2 + Data.Red.Branch.x.^2) .* l;
		tvExpr = (repmat(v, [1 n 1]) - repmat(permute(v, [2 1 3]), [n 1 1])).*Data.Red.Branch.T ...
			- 2 * (Data.Red.Branch.r .* P + Data.Red.Branch.x .* Q) + (Data.Red.Branch.r.^2 + Data.Red.Branch.x.^2) .* l;


		0 >= (tvExpr - repmat(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2, [1 n 1]).*(1-z));
		0 <= (tvExpr + repmat(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2, [1 n 1]).*(1-z));

% 		0 == (1 + repmat(2* Tap .* Data.Red.Bus.Ntr .* Data.Red.Bus.uLow.^2, [1 n 1]) ...
% 			- repmat(permute(v, [2 1 3]), [n 1 1]) + tvEqR).*TSalientesG;

		cDv >= 0;
		cDv >= Data.Cost.m.*(v - (1+Data.Cost.delta));
		cDv >= - Data.Cost.m.*(v - (1-Data.Cost.delta));

		% Restricciones de conexion de red
		Data.Red.Branch.Tup.*z + Data.Red.Branch.Tup.*permute(z, [2 1 3]) == Data.Red.Branch.Tup.*y
		w(G,1,:) == 0;
		0 >= (repmat(w, [1 n]) - repmat(permute(w, [2 1 3]), [n 1 1]) + 1 - 100000*(1-z)).*TnoG;
		w >= 0;
		w <= n-1;


		sum(z(:,NnoG,:).*Data.Red.Branch.T(:,NnoG,:),1) == 1; % los nodos de que no son barras estan todos conectados
		sum(z(:,G,:).*Data.Red.Branch.T(:,G,:),1) == 0; % no hay entradas hacia las barras
		sum(sum(z(G,:,:).*Data.Red.Branch.T(G,:,:))) >= 1; % al menos una barra conectada
		
		Data.Red.Branch.T .* z >= 0;
		Data.Red.Branch.T .* z <= 1;
		
% 		FixY .* y == FixY;
% 		FreeY .* y <= FreeY .* Data.Red.Branch.yTop;
% 		FreeY .* y >= FreeY .* Data.Red.Branch.yLow;
		y <= Data.Red.Branch.yTop;
		y >= Data.Red.Branch.yLow;

		% Restricciones de nodo 0
		% PTras(G,1,:) == sum(Data.Red.Branch.T(G,:,:).*P(G,:,:),2);
		% QTras(G,1,:) == sum(Data.Red.Branch.T(G,:,:).*Q(G,:,:),2);
		% PTras(NnoG,1,:) == 0;
		% QTras(NnoG,1,:) == 0;
		PTras >= Data.Red.Bus.P0Low;
		PTras <= Data.Red.Bus.P0Top;
		QTras >= Data.Red.Bus.Q0Low;
		QTras <= Data.Red.Bus.Q0Top;

		% Restricciones de dominio
		v >= Data.Red.Bus.uLow.^2;
		v <= Data.Red.Bus.uTop.^2;

		
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
		end
			
% 		if Data.Util.Func > 0
        pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnLow.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn; 
        pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnTop.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

        pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnTop.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn;
        pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnLow.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

        qCApp == pCApp.*Data.Util.tgPhi;


        pCn >= Data.Util.pzCnLow;
        pCn <= Data.Util.pzCnTop;

% 		else
% 			pC >= Data.Red.Bus.pCLow.*((1-Data.Red.Bus.alpha(:,1,:,1)) + Data.Red.Bus.alpha.*v);
% 
% 			qC >= Data.Red.Bus.qCLow.*((1-Data.Red.Bus.alpha(:,1,:,1)) + Data.Red.Bus.alpha.*v);
% 
% 		end

		%% Restricciones de generadores Solares
        if lenPv > 0
            pPv == (Data.Gen.Pv.pPvg - (Data.Gen.Pv.cv.*sPv + Data.Gen.Pv.cr.*xiPv)).*Data.Gen.Pv.I;
            pPv == (Data.Gen.Pv.pPvg - (Data.Gen.Pv.cv.*sPv + Data.Gen.Pv.cr.*xiPv));
            sPv >= norms([pPv qPv], 2, 2);
            xiPv >= pPv.^2 + qPv.^2;
            sPv <= Data.Gen.Pv.sTop.*abs(sign(Data.Gen.Pv.pPvg)).*Data.Gen.Pv.I;
            sPv <= Data.Gen.Pv.sTop.*abs(sign(Data.Gen.Pv.pPvg));

            (1-Data.Gen.Pv.I).*pPv == 0;
            (1-Data.Gen.Pv.I).*qPv == 0;
            (1-Data.Gen.Pv.I).*sPv == 0;
            (1-Data.Gen.Pv.I).*xiPv == 0;
        end
        
		%% Restricciones de storage bateria
        if lenSt > 0
            EStbAnt(:,1,1) = Data.St.Bat.EIni(:,1,1);
            EStbAnt(:,1,(2:Config.Etapas)) = EStb(:,1,(1:Config.Etapas-1));

            pStb == pStgb - (Data.St.Bat.cv.*sStb + Data.St.Bat.cr.*xiStb);
            EStb == (1-Data.St.Bat.epsilon).*EStbAnt - Data.St.Bat.eta.*pStgb*Data.dt;
            sStb >= norms([pStgb qStb], 2, 2);
            xiStb >= pStgb.^2 + qStb.^2;

            pStb <= Data.St.Bat.pgTop.*Data.St.Bat.I;
            sStb <= Data.St.Bat.sTop.*Data.St.Bat.I;
            EStb <= Data.St.Bat.ETop.*Data.St.Bat.I;

            pStb >= Data.St.Bat.pgLow.*Data.St.Bat.I;
            EStb >= Data.St.Bat.ELow.*Data.St.Bat.I;
        end

% 		fopt_expr = sum(tfopt_expr) + fopt_expr + tUtilStb;
		fopt_expr = sum(tfopt_expr);

		%% Restricciones de Aire Acondicionado
		TvarAnt(:,1,1) = Data.St.AC.tempIni;
		TvarAnt(:,1,(2:Config.Etapas)) = Tvar(:,1,(1:Config.Etapas-1));
		pCAC(:,:,:) = 0;
		pCAC(iC,1,:) = pCApp(iC,1,:,2);

% 		Tvar(iC,1,:) == TvarAnt(iC,1,:) - Data.St.AC.epsilon(iC,1,:).*(TvarAnt(iC,1,:) - Data.temp(iC,1,:))*Data.dt + Data.St.AC.eta(iC,1,:).*pCAC(iC,1,:)*Data.dt : dTvar;
		Tvar(iC,1,:) - TvarAnt(iC,1,:) + Data.St.AC.epsilon(iC,1,:).*(TvarAnt(iC,1,:) - Data.temp(iC,1,:))*Data.dt - Data.St.AC.eta(iC,1,:).*pCAC(iC,1,:)*Data.dt == 0 : dTvar;
        
		Tvar(iC,1,:) >= Data.St.AC.tempLow(iC,1,:);
		Tvar(iC,1,:) <= Data.St.AC.tempTop(iC,1,:);

		%% Restricciones de generadores Eolico

		% Modelo de Red interna
        if lenWN > 0
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

            vDFIG(2,1,:,:) >= Data.Gen.DFIG.uLow(2,1,:,:).^2;
            vDFIG(2,1,:,:) <= Data.Gen.DFIG.uTop(2,1,:,:).^2;

            vDFIG(3,1,:,:) >= Data.Gen.DFIG.uLow(3,1,:,:).^2;
            vDFIG(3,1,:,:) <= Data.Gen.DFIG.uTop(3,1,:,:).^2;

            Data.Gen.DFIG.PQnorm(1,2,:,:) >= norms([PDFIG(1,2,:,:) QDFIG(1,2,:,:)],2 ,2);
            Data.Gen.DFIG.PQnorm(1,3,:,:) >= norms([PDFIG(1,3,:,:) QDFIG(1,3,:,:)],2 ,2);

            sDFIG (3,1,:,:) >= norms([pCDFIG(3,1,:,:) qCDFIG(3,1,:,:)],2 ,2);
            xiDFIG(3,1,:,:) >= pCDFIG(3,1,:,:).^2 + qCDFIG(3,1,:,:).^2;

            sDFIG (5,1,:,:) >= norms([PDFIG(4,5,:,:) QDFIG(4,5,:,:)],2 ,2);
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

        end


        %% Restricciones de cargas no interrumpibles
        if lenCh > 0
            onNext(nodCh,1,(1:Config.Etapas-1)) = onCh(nodCh,1,(2:Config.Etapas));
            onNext(nodCh,1,Config.Etapas) = 0;
            stCh - onNext + onCh >= 0;

            sum(stCh(nodCh,1,:),3) == 1;
            sum(onCh(nodCh,1,:),3) == Data.ClNI.d(nodCh);
            stCh(nodCh,1,1) == onCh(nodCh,1,1);

            stCh(NnodCh,1,:) == 0;
            onCh(NnodCh,1,:) == 0;
        end
        
        if fixed
            Cap == Data.Fixed.Cap;
            Tap == Data.Fixed.Tap;
            y == Data.Fixed.y;
        end
        if lenCh > 0
            if fixed
                stCh == Data.Fixed.stCh;
                onCh == Data.Fixed.onCh;
            end
        end

		if findOpt
			minimize fopt_expr
        end
        
    cvx_end

	toc
	
    if lenCh > 0
        Var.ClNI.pC = Data.ClNI.pC.*onCh;
        Var.ClNI.qC = Data.ClNI.qC.*onCh;
        Var.ClNI.on = onCh;
        Var.ClNI.start = stCh;
        Var.Cost.optEt = tfopt_expr;
    end

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
	else
		Var.Gen.Dfig.pWi = zeros(n,Config.Etapas);
		Var.Gen.Dfig.qWi = zeros(n,Config.Etapas);
	end

    Var.Red.Branch.P = P;
    Var.Red.Branch.Q = Q;
	Var.Red.Branch.l = l;
	Var.Red.Branch.lNorm = lNorm;
	Var.Red.Branch.lQoL = lQoL;
    Var.Red.Bus.w = w;
	Var.Red.Branch.y = y;
	Var.Red.Branch.z = z;

%     for i = 1:size(stCh,1)
% 		dur = ones(1,Cargas(i).dur);
% 		ini = find(Var.ClNI.on(i,:) == 1);
% 		Var.ClNI.on(i,ini:ini+Cargas(i).dur-1) = dur;
% 	end
	Var.Red.Bus.Cap = Cap;
	Var.Red.Bus.PTras = PTras;
	Var.Red.Bus.QTras = QTras;
	Var.Red.Bus.Tap = Tap;
	Var.Red.Bus.cDv = cDv;
	Var.Red.Bus.nn = nn;
	Var.Red.Bus.nv = nv;
% 	Var.Red.Bus.pC = pC;
% 	Var.Red.Bus.pCApp = pCApp;
	Var.Red.Bus.pN = pN;
% 	Var.Red.Bus.qC = qC;
% 	Var.Red.Bus.qCApp = qCApp;
	Var.Red.Bus.qCp = qCp;
	Var.Red.Bus.qN = qN;
	Var.Red.Bus.v = v;

	Var.ClRes.pC = sum(pCApp,4);
	Var.ClRes.qC = sum(qCApp,4);
	Var.ClRes.pCApp = pCApp;
	Var.ClRes.qCApp = qCApp;
	Var.ClRes.Tvar = Tvar;

    
	
	% auxCapLow = repmat(Data.Red.Bus.CapLow,1,Config.Etapas);
	% auxCapTop = repmat(Data.Red.Bus.CapTop,1,Config.Etapas);
	% auxULow = repmat(Data.Red.Bus.uLow,1,Config.Etapas);
	% auxUTop = repmat(Data.Red.Bus.uTop,1,Config.Etapas);

%     Var.Red.Bus.mcCaps = max(auxCapLow.*v + Cap.*(auxULow).^2 - auxCapLow.*(auxULow).^2, ...
%         auxCapTop.*v + Cap.*(auxUTop).^2 - auxCapTop.*(auxUTop).^2);

    
	% if Data.Util.Func > 0
		% Var.Red.Bus.pzC = pCApp;
		% Var.Red.Bus.pCn = pCn;
	% else
		% Var.Red.Bus.pzC = zeros(n, Config.Etapas, 2);
		% Var.Red.Bus.pCn = zeros(n, Config.Etapas, 2);
	% end

    if fixed
        Var.Dual.ind = NnoG;
        Var.Dual.dPn = dPn;
        Var.Dual.dQn = dQn;
        Var.Dual.dTvar = dTvar;
%		 Var.Dual.dQnT = cell2mat(dQnT');
    end
    
    if lenPv > 0
        Var.Gen.Pv.CqgS = sum(sum(cqPv,1),3);
        Var.Gen.Pv.pPv = pPv;
        Var.Gen.Pv.qPv = qPv;
        Var.Gen.Pv.s = sPv;
        Var.Gen.Pv.xi = xiPv;
    end

% 	Var.Cost.CpgS = Data.Cost.rhopPv' * pPv;
% 	Var.Cost.CpgW = Data.Cost.rhopWi' * pWi;
	% Var.Cost.CqgW = sum(cqWi);
% 	Var.Cost.CPn = Data.Cost.piPTras' * PTras;
	% Var.Cost.CQn = sum(cQTras);
% 	Var.Cost.CCDv = Data.Cost.cdv' * cDv;
	% Var.Cost.CUtil = - sum(tUtil);

% 	Var.Opt.CTras = sum(tCTras);
% 	Var.Opt.CDV = Data.Cost.cdv' * cDv;

if lenSt > 0
	Var.St.AC.Tvar = Tvar;
	Var.St.Bat.pStb = pStb;
	Var.St.Bat.pStgb = pStgb;
	Var.St.Bat.qStb = qStb;
	Var.St.Bat.sStb = sStb;
	Var.St.Bat.xiStb = xiStb;
	Var.St.Bat.EStb = EStb;
end
    
    
	redCost = 0;
	if isfield(Data.Red, 'Cost')
		redCost = Data.Red.Cost;
	end
	opt = cvx_optval + redCost;
    
    Var.Red.Bus.pG = pG;
	Var.Red.Bus.qG = qG;

	Var.Red.Bus.pC = pC;
	Var.Red.Bus.qC = qC;

	status = cvx_status;

	% if (lenWN > 0)
		% Var.Gen.Dfig.Branch.errLR = Var.Gen.Dfig.Branch.l * 0;
		% Var.Gen.Dfig.Branch.errLA = Var.Gen.Dfig.Branch.errLR;

%         for i = 1:Config.Etapas
%             PQV = (Var.Gen.Dfig.Branch.P(:,:,i).^2 + Var.Gen.Dfig.Branch.Q(:,:,i).^2) ./ (Var.Gen.Dfig.Bus.v(:,:,i)'*ones(1,length(Data.Gen.DFIG.Tg)).*Data.Gen.DFIG.Tg);
%             Var.Gen.Dfig.Branch.errLR(:,:,i) = abs(Var.Gen.Dfig.Branch.l(:,:,i) - PQV)./PQV;
%             Var.Gen.Dfig.Branch.errLA(:,:,i) = Data.Gen.DFIG.r .* abs(Var.Gen.Dfig.Branch.l(:,:,i) - PQV);
%         end

		% Var.Gen.Dfig.Bus.errSR = Var.Gen.Dfig.Bus.s * 0;
		% Var.Gen.Dfig.Bus.errSA = Var.Gen.Dfig.Bus.errSR;
		% Var.Gen.Dfig.Bus.errXiR = Var.Gen.Dfig.Bus.errSR;
		% Var.Gen.Dfig.Bus.errXiA = Var.Gen.Dfig.Bus.errSR;
		% auxSDFIG = Var.Gen.Dfig.Bus.errSR;
		% auxSDFIG(3,1,:) = Var.Gen.Dfig.Bus.pC(3,1,:).^2 + Var.Gen.Dfig.Bus.qC(3,1,:).^2;
		% auxSDFIG(5,1,:) = Var.Gen.Dfig.Branch.P(4,5,:).^2 + Var.Gen.Dfig.Branch.Q(4,5,:).^2;

		% Var.Gen.Dfig.Bus.errSR = abs(Var.Gen.Dfig.Bus.s - sqrt(auxSDFIG)) ./ sqrt(auxSDFIG);
		% Var.Gen.Dfig.Bus.errXiR = abs(Var.Gen.Dfig.Bus.xi - auxSDFIG) ./ auxSDFIG;

%         for i = 1:Config.Etapas
%             Var.Gen.Dfig.Bus.errSA(:,:,i) = Data.Gen.DFIG.cv' .* abs(Var.Gen.Dfig.Bus.s(:,:,i) - sqrt(auxSDFIG(:,:,i)));
% 
%             Var.Gen.Dfig.Bus.errXiA(:,:,i) = Data.Gen.DFIG.cr' .* abs(Var.Gen.Dfig.Bus.xi(:,:,i) - auxSDFIG(:,:,i));
%         end
	% end

	% for et = 1:Config.Etapas
    
%         for app = 1:2
%             LL = Data.Util.pzCnLow(:,et,app).*Var.Red.Bus.v(:,app) + pCn(:,1,et,app).*Data.Red.Bus.uLow.^2 - Data.Util.pzCnLow(:,et,app).*Data.Red.Bus.uLow.^2;
%             TT = Data.Util.pzCnTop(:,et,app).*Var.Red.Bus.v(:,app) + pCn(:,1,et,app).*Data.Red.Bus.uTop.^2 - Data.Util.pzCnTop(:,et,app).*Data.Red.Bus.uTop.^2;
%             TL = Data.Util.pzCnTop(:,et,app).*Var.Red.Bus.v(:,app) + pCn(:,1,et,app).*Data.Red.Bus.uLow.^2 - Data.Util.pzCnTop(:,et,app).*Data.Red.Bus.uLow.^2;
%             LT = Data.Util.pzCnLow(:,et,app).*Var.Red.Bus.v(:,app) + pCn(:,1,et,app).*Data.Red.Bus.uTop.^2 - Data.Util.pzCnLow(:,et,app).*Data.Red.Bus.uTop.^2;
% 
%             LL = Data.Red.Bus.alpha(:,app).*(LL) + (1-Data.Red.Bus.alpha(:,app)).* pCn(:,1,et,app);
%             TT = Data.Red.Bus.alpha(:,app).*(TT) + (1-Data.Red.Bus.alpha(:,app)).* pCn(:,1,et,app);
%             TL = Data.Red.Bus.alpha(:,app).*(TL) + (1-Data.Red.Bus.alpha(:,app)).* pCn(:,1,et,app);
%             LT = Data.Red.Bus.alpha(:,app).*(LT) + (1-Data.Red.Bus.alpha(:,app)).* pCn(:,1,et,app);
% 
%             Var.Red.Bus.errLL(:,et,app) = abs((pCApp(:,1,et,app) - LL) ./ LL);
%             Var.Red.Bus.errTT(:,et,app) = abs((pCApp(:,1,et,app) - TT) ./ TT);
%             Var.Red.Bus.errTL(:,et,app) = abs((pCApp(:,1,et,app) - TL) ./ TL);
%             Var.Red.Bus.errLT(:,et,app) = abs((pCApp(:,1,et,app) - LT) ./ LT);
%         end
		% Var.Red.Bus.elastP(:,et) = abs((pCn(:,1,et,1) - Data.Util.pzCnPref(:,1)) ./ Data.Util.pzCnPref(:,1));
        % Var.Red.Bus.elastT(:,et) = abs((Tvar(:,1,et) - Data.St.AC.tempPref) ./ Data.St.AC.tempPref);
	% end

	% auxSPv = Var.Gen.Pv.pPv.^2 + Var.Gen.Pv.qPv.^2;
	% Var.Gen.Pv.errS = abs(Var.Gen.Pv.s - sqrt(auxSPv)) ./ sqrt(auxSPv);
	% Var.Gen.Pv.errXi = abs(Var.Gen.Pv.s - auxSPv) ./ auxSPv;

	% Var.Red.Branch.errLR = Var.Red.Branch.l * 0;
	% Var.Red.Branch.errLA = Var.Red.Branch.errLR;

%     for i = 1:Config.Etapas
%         PQV = (Var.Red.Branch.P(:,:,i).^2 + Var.Red.Branch.Q(:,:,i).^2) ./ ((Var.Red.Bus.v(:,1,i)*ones(1,length(Data.Red.Branch.T(:,:,i)).*Data.Red.Branch.T(:,:,i)).*Var.Red.Branch.z(:,:,i)));
%         Var.Red.Branch.errLR(:,:,i) = Var.Red.Branch.z(:,:,i).*abs(Var.Red.Branch.l(:,:,i) - PQV)./PQV;
%         Var.Red.Branch.errLA(:,:,i) = Var.Red.Branch.z(:,:,i).*Data.Red.Branch.r.*abs(Var.Red.Branch.l(:,:,i) - PQV);
%     end
    
    
    cvx_clear;
	end;
end
