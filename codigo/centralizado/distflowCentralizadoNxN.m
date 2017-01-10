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




NoT = 1 - Data.Red.Branch.T;

tnnLow = (1 + Data.Red.Bus.TapLow.*Data.Red.Bus.Ntr);
tnnTop = (1 + Data.Red.Bus.TapTop.*Data.Red.Bus.Ntr);

%% Inicializacion

nodPv = find(Data.Gen.Pv.I == 1);
lenPv = length(nodPv);

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
        
        variable pC(n,1, Config.Etapas); % real power demand in i
        variable qC(n,1, Config.Etapas); % reactive power demand in i
        
        variable pN(n,1, Config.Etapas);
        variable qN(n,1, Config.Etapas);
        
        variable pG(n,1, Config.Etapas);
        variable qG(n,1, Config.Etapas);
		
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

		variable Cap(n,1, Config.Etapas) integer;
		variable Tap(n,1, Config.Etapas) integer;
		variable y(n, n, Config.Etapas) binary;
        
        
		% Variables de utilidad
        variable pCApp(n,1, Config.Etapas, 2); % real power demand in i
        variable qCApp(n,1, Config.Etapas, 2); % real power demand in i
        expression vApp(n, 1, Config.Etapas, 2);
        variable pCn(n,1, Config.Etapas, 2); % real power demand in i



        variable pCClRes(n,1, Config.Etapas); % real power demand in i
        variable qCClRes(n,1, Config.Etapas); % real power demand in i
      
        

		%% Declaracion de variables y expresiones
		expression tvExpr(n,n,Config.Etapas);
		expression CapDif(n,1,Config.Etapas);
		expression TapDif(n,1,Config.Etapas);

		% Variables de generadores solares
        variable pPv(n, 1, Config.Etapas); % real power generation in i
        variable qPv(n, 1, Config.Etapas); % reactive power generation in i


		expression tfopt_expr(Config.Etapas); 
		expression fopt_expr; 

        
		expression lQoL(n, n, Config.Etapas,3);
		expression lNorm(n, n, Config.Etapas,3);


		%% Funcion objetivo
        CapDif(:,1,1) = Data.Red.Bus.CapIni;
        CapDif(:,1,(2:Config.Etapas)) = Cap(:,1,(2:Config.Etapas)) - Cap(:,1,(1:Config.Etapas-1));

        TapDif(:,1,1) = Data.Red.Bus.TapIni;
        TapDif(:,1,(2:Config.Etapas)) = Tap(:,1,(2:Config.Etapas)) - Tap(:,1,(1:Config.Etapas-1));
		
		tfopt_expr = ...
            sum(Data.Cost.piPTras .* PTras,1) ...
            + sum(Data.Cost.cdv .* cDv,1) ...
            + sum(cQTras,1) ...
            + sum(Data.Red.cambioCap*CapDif.^2,1) ...
            + sum(Data.Red.cambioTap*TapDif.^2,1) ...
            + sum(Data.Util.betaT(:,1,:,1).*(pCn(:,1,:,1) - Data.Util.pzCnPref(:,1,:,1)).^2,1) ...
            ;


        
        
        

		cQTras >= - Data.Cost.piQmtras .* QTras;
		cQTras >= Data.Cost.piQMtras .* QTras;


		%% Restricciones de Red
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


        pC == pCClRes;
        qC == qCClRes;
        
        pG == PTras + pPv;
        qG == QTras + qPv + qCp;
        
		pN - pC + pG == 0;
		qN - qC + qG == 0;

		% Restricciones de la corriente
        lQoL(:,:,:,1) = Data.Red.Branch.T.*(2*P);
        lQoL(:,:,:,2) = Data.Red.Branch.T.*(2*Q);
        lQoL(:,:,:,3) =  Data.Red.Branch.T.*(l - repmat(v, [1 n 1]));
        norms(lQoL,2,4) <= Data.Red.Branch.T.*(l + repmat(v, [1 n 1]));

        lNorm(:,:,:,1) = Data.Red.Branch.T.*(2*P);
        lNorm(:,:,:,2) = Data.Red.Branch.T.*(2*Q);
        lNorm(:,:,:,3) =  Data.Red.Branch.T.*(l - repmat(Data.Red.Bus.uTop.^2, [1 n 1]).*z);
        norms(lNorm,2,4) <= Data.Red.Branch.T.*(l + repmat(Data.Red.Bus.uTop.^2, [1 n 1]).*z);

		l <= Data.Red.Branch.lTop.*z;

		% Restricciones de voltaje

		nn >= (1 + Tap.*Data.Red.Bus.Ntr).^2;
		nn <= (tnnTop + tnnLow).*(1 + Tap.*Data.Red.Bus.Ntr) - (tnnTop.*tnnLow);

		nv >= nn.*(Data.Red.Bus.uLow.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uLow.^2);
		nv >= nn.*(Data.Red.Bus.uTop.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uTop.^2);

		nv <= nn.*(Data.Red.Bus.uLow.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uLow.^2);
		nv <= nn.*(Data.Red.Bus.uTop.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uTop.^2);

		Tap >= Data.Red.Bus.TapLow;
		Tap <= Data.Red.Bus.TapTop;

		tvExpr = (repmat(nv, [1 n 1]) - repmat(permute(v, [2 1 3]), [n 1 1])).*Data.Red.Branch.T ...
			- 2 * (Data.Red.Branch.r .* P + Data.Red.Branch.x .* Q) + (Data.Red.Branch.r.^2 + Data.Red.Branch.x.^2) .* l;


        0 >= (tvExpr - Data.Red.Branch.T.*repmat(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2, [1 n 1]).*(1-z));
        0 <= (tvExpr + Data.Red.Branch.T.*repmat(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2, [1 n 1]).*(1-z));

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
		
		y <= Data.Red.Branch.yTop;
		y >= Data.Red.Branch.yLow;

		% Restricciones de nodo 0
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
			
        pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnLow.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn; 
        pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnTop.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

        pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnTop.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn;
        pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnLow.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

        qCApp == pCApp.*Data.Util.tgPhi;


        pCn >= Data.Util.pzCnLow;
        pCn <= Data.Util.pzCnTop;

        pCClRes >= sum(pCApp, 4);

        qCClRes >= sum(qCApp, 4);

		%% Restricciones de generadores Solares
		% Variables de generadores solares
        if lenPv > 0
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
            
            pPv == (Data.Gen.Pv.pPvg - (Data.Gen.Pv.cv.*sPv + Data.Gen.Pv.cr.*xiPv)).*Data.Gen.Pv.I;

            SPvNorm(:,:,:,1) = pPv;
            SPvNorm(:,:,:,2) = qPv;
            sPv >= norms(SPvNorm,2,4);
            xiPv >= pPv.^2 + qPv.^2;
            
            sPv <= Data.Gen.Pv.sTop.*abs(sign(Data.Gen.Pv.pPvg)).*Data.Gen.Pv.I;

        else
            pPv == 0;
            qPv == 0;
        end
        
		fopt_expr = sum(tfopt_expr);

		if findOpt
			minimize fopt_expr
        end
        
    cvx_end

	toc
	

    Var.Red.Branch.P = P;
    Var.Red.Branch.Q = Q;
	Var.Red.Branch.l = l;
	Var.Red.Branch.lNorm = lNorm;
	Var.Red.Branch.lQoL = lQoL;
    Var.Red.Bus.w = w;
	Var.Red.Branch.y = y;
	Var.Red.Branch.z = z;

	Var.Red.Bus.Cap = Cap;
	Var.Red.Bus.PTras = PTras;
	Var.Red.Bus.QTras = QTras;
	Var.Red.Bus.Tap = Tap;
	Var.Red.Bus.cDv = cDv;
	Var.Red.Bus.nn = nn;
	Var.Red.Bus.nv = nv;
	Var.Red.Bus.pN = pN;
	Var.Red.Bus.qCp = qCp;
	Var.Red.Bus.qN = qN;
	Var.Red.Bus.v = v;

    Var.ClRes.pC = pCClRes;
    Var.ClRes.qC = qCClRes;
	Var.ClRes.pCApp = pCApp;
	Var.ClRes.qCApp = qCApp;

    if lenPv > 0
        Var.Gen.Pv.CqgS = sum(sum(cqPv,1),3);
        Var.Gen.Pv.pPv = pPv;
        Var.Gen.Pv.qPv = qPv;
        Var.Gen.Pv.s = sPv;
        Var.Gen.Pv.xi = xiPv;
    end

	opt = cvx_optval;
    
    Var.Red.Bus.pG = pG;
	Var.Red.Bus.qG = qG;

	Var.Red.Bus.pC = pC;
	Var.Red.Bus.qC = qC;

	status = cvx_status;

    
    cvx_clear;
end
