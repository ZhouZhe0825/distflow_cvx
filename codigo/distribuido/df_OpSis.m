function [Var, opt, status] = df_OpSis(Data, Config, findOpt, DistrInfo)
% function [Var, opt, status] = distflow_nod_eol_util(Data, rhopgSn, rhoP, rhopgWn, rhomqgSn, rhoMqgSn, rhomQ, rhoMQ, rhomqgWn, rhoMqgWn, delta, m, cdv, findOpt, fijarPgSn)
% Data.Red.Branch.T, Data.Red.Branch.r, Data.Red.Branch.x, Data.Gen.Pv.I son matrices cuadradas n * n ...
% Data.Red.Bus.alpha, Data.Red.Bus.cv, Data.Red.Bus.cr, Data.Gen.Pv.qgTop, uTop, Data.Red.Bus.uLow, Data.Red.Bus.pcLow, Data.Red.Bus.qcLow, Data.Red.Bus.qcCap, Data.Gen.Pv.sTop, Data.Gen.Pv.pgTop son vectores columnas de largo n
% Data.Red.Bus.v0, Data.Red.Bus.Q0Top, Data.Red.Bus.Q0Low, Data.Red.Bus.P0Top, Data.Red.Bus.P0Low son escalares
% modSym, findOpt son booleano

    fixed = isfield(Data, 'Fixed');



	n = size(Data.Red.Branch.T,1);

	%% Inicializacion
	G = Data.Red.Bus.v0;
	NnoG = setdiff((1:n), G)';

	TSalientesG = Data.Red.Branch.T;
	TEntrantesG = Data.Red.Branch.T;
	TnoG = Data.Red.Branch.T;
	TnoG(G,:,:) = zeros(length(G),n,Config.Etapas);
	TSalientesG = TSalientesG - TnoG;
	TnoG(:,G,:) = zeros(n,length(G),Config.Etapas);
	TEntrantesG = TEntrantesG - TSalientesG - TnoG;

	FixY = Data.Red.Branch.T .* Data.Red.Branch.yTop .* Data.Red.Branch.yLow;
	FreeY = Data.Red.Branch.T .* (1 - (Data.Red.Branch.yTop .* Data.Red.Branch.yLow));
	NoT = 1 - Data.Red.Branch.T;

	NcpCapL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow;
	NcpCapT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop;

	NcpvL = Data.Red.Bus.Ncp.*(Data.Red.Bus.uLow).^2;
	NcpvT = Data.Red.Bus.Ncp.*(Data.Red.Bus.uTop).^2;

	NcpCapLvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uLow).^2;
	NcpCapTvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uTop).^2;
	NcpCapLvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uTop).^2;
	NcpCapTvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uLow).^2;

	tnnLow = (1 + Data.Red.Bus.TapLow.*Data.Red.Bus.Ntr);
	tnnTop = (1 + Data.Red.Bus.TapTop.*Data.Red.Bus.Ntr);

	isSymT = isequal(matOverTime(Data.Red.Branch.T), matOverTime(Data.Red.Branch.T)');

	if isSymT
		%% Modelo programacion matematica

		cvx_begin quiet

		% solo para mosek

            for i = 1: size(Config,1)
                cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
            end

			%% Declaracion de variables y expresiones
			%Variables generales de la red
			variable P(n, n, Config.Etapas); % active power from i to j
			variable Q(n, n, Config.Etapas); % reactive power from i to j
			variable l(n, n, Config.Etapas);
			variable pN(n, 1, Config.Etapas); % real power demand in i
			variable qN(n, 1, Config.Etapas); % reactive power demand in i
			variable qCp(n, 1, Config.Etapas); % reactive power demand in i
			variable v(n, 1, Config.Etapas); % module of square complex voltage in i
			variable w(n, 1, Config.Etapas);
			variable nn(n, 1, Config.Etapas);
			variable nv(n, 1, Config.Etapas);

			variable cDv(n, 1, Config.Etapas);

%             variable z(n, n, Config.Etapas);
            if fixed
                variable Cap(n, 1, Config.Etapas);
                variable Tap(n, 1, Config.Etapas);
                variable y(n, n, Config.Etapas);
            else
                variable Cap(n, 1, Config.Etapas) integer;
                variable Tap(n, 1, Config.Etapas) integer;
                variable y(n, n, Config.Etapas) binary;
            end
            variable z(n, n, Config.Etapas);
            
			% Expresion general de (P, pg's, pc) y (Q, qg's, qc) de la red
			expression tvEqL(n,n,Config.Etapas);
			expression tvEqR(n,n,Config.Etapas);
			expression CapDif(n,1,Config.Etapas);
			expression TapDif(n,1,Config.Etapas);

			% Expresion para p y q generados que entran en la red general
			expression tPT(n, n, Config.Etapas);
			expression tQT(n, n, Config.Etapas);

			expression tfopt_expr(Config.Etapas);
			expression tfopt_virt(Config.Etapas);
			expression tfopt_conv(Config.Etapas);
			expression fopt_expr; 

			expression lQoL(n, n, Config.Etapas,3);
			expression lNorm(n, n, Config.Etapas,3);

			fopt_expr = 0;

			%% Modelado de etapas

			tPT = Data.Red.Branch.T.*P;
			tQT = Data.Red.Branch.T.*Q;


			%% Restricciones de Red

			% Restricciones generales de branches
			pN == (permute(sum(tPT - Data.Red.Branch.r.*l, 1),[2 1 3]) - sum(tPT, 2));
			qN == (permute(sum(tQT - Data.Red.Branch.x.*l, 1),[2 1 3]) - sum(tQT, 2));
            
            pN(Data.Red.Bus.indCG0,:,:) == 0;
            qN(Data.Red.Bus.indCG0,:,:) == 0;
            

			Cap >= Data.Red.Bus.CapLow;
			Cap <= Data.Red.Bus.CapTop;

			% Restricciones de potencias de la red
% % % % 			if fixed
% % % % 				qCp == Data.Red.Bus.Ncp.*Data.Fixed.Cap.*v;
% % % % 			else
                qCp >= NcpCapL.*v + NcpvL.*Cap - NcpCapLvL;
                qCp >= NcpCapT.*v + NcpvT.*Cap - NcpCapTvT;
                qCp <= NcpCapL.*v + NcpvT.*Cap - NcpCapLvT;
                qCp <= NcpCapT.*v + NcpvL.*Cap - NcpCapTvL;
% % % % 			end

			%% Funcion objetivo
            CapDif(:,1,1) = Data.Red.Bus.CapIni;
            CapDif(:,1,(2:Config.Etapas)) = Cap(:,1,(2:Config.Etapas)) - Cap(:,1,(1:Config.Etapas-1));

            TapDif(:,1,1) = Data.Red.Bus.TapIni;
            TapDif(:,1,(2:Config.Etapas)) = Tap(:,1,(2:Config.Etapas)) - Tap(:,1,(1:Config.Etapas-1));

            tfopt_expr = sum(Data.Cost.cdv .* cDv + Data.Red.cambioCap*CapDif(:,1,:).^2 + Data.Red.cambioTap*TapDif(:,1,:).^2,1);
            tfopt_virt = sum(- DistrInfo.muT(:,1,:) .* pN - DistrInfo.lambdaT(:,1,:) .* qN + DistrInfo.lambdaT(:,1,:) .* qCp,1);
            tfopt_conv = sum(norms(pN - DistrInfo.OpSis.pN,2,2) ...
                            + norms(qN - DistrInfo.OpSis.qN,2,2) ...
                            + norms(qCp - DistrInfo.OpSis.qCp,2,2),1);
%             if fixed
                tfopt_conv = 1/(2*DistrInfo.Gama) * tfopt_conv;
%                 tfopt_conv = 25000 * tfopt_conv;
%             else
%                 tfopt_conv = 1/(2*DistrInfo.Gama) * tfopt_conv;
% %                 tfopt_conv = 50000 * tfopt_conv;
%             end
                fopt_expr = sum(tfopt_expr) + sum(tfopt_virt) + sum(tfopt_conv);

			% Restricciones de la corriente
            lQoL(:,:,:,1) = 2*P.*Data.Red.Branch.T;
            lQoL(:,:,:,2) = 2*Q.*Data.Red.Branch.T;
            lQoL(:,:,:,3) = (l - repmat(v, [1 n 1])).*Data.Red.Branch.T;
            norms(lQoL,2,4) - (l + repmat(v, [1 n 1])).*Data.Red.Branch.T <= 0;

            lNorms(:,:,:,1) = 2*P.*Data.Red.Branch.T;
            lNorms(:,:,:,2) = 2*Q.*Data.Red.Branch.T;
            lNorms(:,:,:,3) = (l- z.*repmat(Data.Red.Bus.uTop.^2, [1 n 1])).*Data.Red.Branch.T;
            norms(lNorms,2,4) - (l + z.*repmat(Data.Red.Bus.uTop.^2, [1 n 1])).*Data.Red.Branch.T <= 0;
 			
			% Restricciones de voltaje
			nn >= (1 + Tap.*Data.Red.Bus.Ntr(:,1,:)).^2;
			nn <= (tnnTop + tnnLow).*(1 + Tap.*Data.Red.Bus.Ntr(:,1,:)) - (tnnTop.*tnnLow);

			nv >= nn.*(Data.Red.Bus.uLow.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uLow.^2);
			nv >= nn.*(Data.Red.Bus.uTop.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uTop.^2);

			nv <= nn.*(Data.Red.Bus.uLow.^2) + tnnTop.*v - tnnTop.*(Data.Red.Bus.uLow.^2);
			nv <= nn.*(Data.Red.Bus.uTop.^2) + tnnLow.*v - tnnLow.*(Data.Red.Bus.uTop.^2);

			Tap >= Data.Red.Bus.TapLow(:,1,:);
			Tap <= Data.Red.Bus.TapTop(:,1,:);

			tvEqR = - 2 * (Data.Red.Branch.r .* P + Data.Red.Branch.x .* Q) + (Data.Red.Branch.r.^2 + Data.Red.Branch.x.^2) .* l;
			tvEqL = (repmat(nv, [1 n 1]) - repmat(permute(v, [2 1 3]), [n 1 1])).*Data.Red.Branch.T;
			0 >= (tvEqL + tvEqR - repmat(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2, [1 n 1]).*(1-z)).*TnoG;
			0 <= (tvEqL + tvEqR + repmat(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2, [1 n 1]).*(1-z)).*TnoG;
% 			0 == (repmat(nn.*Data.Red.Bus.uLow.^2, [1 n 1]) ...
% 				- repmat(permute(v, [2 1 3]), [n 1 1]) + tvEqR).*TSalientesG;
			
			0 == (1 + repmat(2* Tap .* Data.Red.Bus.Ntr .* Data.Red.Bus.uLow.^2, [1 n 1]) ...
				- repmat(permute(v, [2 1 3]), [n 1 1]) + tvEqR).*TSalientesG;

            cDv >= 0;
			cDv >= Data.Cost.m.*(v - (1+Data.Cost.delta));
			cDv >= - Data.Cost.m.*(v - (1-Data.Cost.delta));

			% Restricciones de conexion de red
            Data.Red.Branch.Tup.*z + Data.Red.Branch.Tup.*permute(z, [2 1 3]) == Data.Red.Branch.Tup.*y
            w(G,1,:) == 0;
            w >= 0;
            w <= n-1;

            sum(z(:,NnoG,:).*Data.Red.Branch.T(:,NnoG,:),1) == 1; % los nodos de que no son barras estan todos conectados
            sum(z(:,G,:).*Data.Red.Branch.T(:,G,:),1) == 0; % no hay entradas hacia las barras
            sum(sum(z(G,:,:).*Data.Red.Branch.T(G,:,:))) >= 1; % al menos una barra conectada

            Data.Red.Branch.T .* z >= 0;
            Data.Red.Branch.T .* z <= 1;

            FixY .* y == FixY;
            FreeY .* y <= FreeY .* Data.Red.Branch.yTop;
            FreeY .* y >= FreeY .* Data.Red.Branch.yLow;

			% Restricciones de dominio
			v >= (Data.Red.Bus.uLow).^2;
			v <= (Data.Red.Bus.uTop).^2;

			Data.Red.Branch.T.*l <= Data.Red.Branch.T.*Data.Red.Branch.lTop.*z;

			NoT.*P == 0;
			NoT.*Q == 0;
			NoT.*l == 0;
			NoT.*z == 0;
			NoT.*y == 0;

			P(Data.Red.Bus.v0,:,:) == DistrInfo.OpSis.P(Data.Red.Bus.v0,:,:);
			Q(Data.Red.Bus.v0,:,:) == DistrInfo.OpSis.Q(Data.Red.Bus.v0,:,:);


            if fixed
                Cap == Data.Fixed.Cap;
                Tap == Data.Fixed.Tap; 
                y == Data.Fixed.y;
                z == Data.Fixed.z;
            end
            
            
			if findOpt
				minimize fopt_expr
			end
			

		cvx_end

		Var.Red.Branch.P = P;
		Var.Red.Branch.Q = Q;
		Var.Red.Branch.z = z;
		Var.Red.Branch.y = y;
		Var.Red.Bus.w = w;

		Var.Red.Bus.pN = pN;
		Var.Red.Bus.qN = qN;
		Var.Red.Bus.qCp = qCp;
		Var.Red.Bus.v = v;
		Var.Red.Bus.cDv = cDv;
		Var.Red.Branch.l = l;
		Var.Red.Bus.Cap = Cap;
		Var.Red.Bus.Tap = Tap;

		Var.Cost.optEt = tfopt_expr;

		Var.Red.Branch.lQoL = lQoL;
		Var.Red.Branch.lNorm = lNorm;
		
		Var.Red.Bus.nn = nn;
		Var.Red.Bus.nv = nv;
		
		redCost = 0;
		if isfield(Data.Red, 'Cost')
			redCost = Data.Red.Cost;
        end

		% Var.errCapPrin = abs(qCp - Cap.*Data.Red.Bus.Ncp.*v) ./ (abs(Cap.*Data.Red.Bus.Ncp.*v) + eps);
		% Var.errCapLL = abs(qCp - NcpCapL.*v + NcpvL.*Cap - NcpCapLvL) ./ (abs(NcpCapL.*v + NcpvL.*Cap - NcpCapLvL) + eps);
		% Var.errCapTT = abs(qCp - NcpCapT.*v + NcpvT.*Cap - NcpCapTvT) ./ (abs(NcpCapT.*v + NcpvT.*Cap - NcpCapTvT) + eps);
		% Var.errCapLT = abs(qCp - NcpCapL.*v + NcpvT.*Cap - NcpCapLvT) ./ (abs(NcpCapL.*v + NcpvT.*Cap - NcpCapLvT) + eps);
		% Var.errCapTL = abs(qCp - NcpCapT.*v + NcpvL.*Cap - NcpCapTvL) ./ (abs(NcpCapT.*v + NcpvL.*Cap - NcpCapTvL) + eps);

		Var.Red.Bus.qG = Var.Red.Bus.qCp;
        
        opt(1,1) = sum(tfopt_expr) + redCost;
		opt(1,2) = sum(tfopt_virt);
		opt(1,3) = sum(tfopt_conv);
		status = cvx_status;

		% Var.Red.Branch.errLR = Var.Red.Branch.l * 0;
		% Var.Red.Branch.errLA = Var.Red.Branch.errLR;

		% PQV = (Var.Red.Branch.P.^2 + Var.Red.Branch.Q.^2) ./ ((repmat(Var.Red.Bus.v, [1 n 1]).*Data.Red.Branch.T).*Var.Red.Branch.z);
		% Var.Red.Branch.errLR = Var.Red.Branch.z.*abs(Var.Red.Branch.l - PQV)./PQV;
		% Var.Red.Branch.errLA = Var.Red.Branch.z.*Data.Red.Branch.r.*abs(Var.Red.Branch.l - PQV);

		% Var.Red.Bus.errn = abs(nn - (1 + Tap.*Data.Red.Bus.Ntr(:,1,:)).^2) ./ (1 + Tap.*Data.Red.Bus.Ntr(:,1,:));

        
        cvx_clear;
	end;
end
