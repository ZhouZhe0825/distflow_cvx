function [Var, opt, status] = df_ClRes(Data, Config, findOpt, DistrInfo)
% function [Var, opt, status] = distflow_nod_eol_util(Data, rhopgSn, rhoP, rhopgWn, rhomqgSn, rhoMqgSn, rhomQ, rhoMQ, rhomqgWn, rhoMqgWn, delta, m, cdv, findOpt, fijarPgSn)
% Data.Red.Branch.T, Data.Red.Branch.r, Data.Red.Branch.x, Data.Gen.Pv.I son matrices cuadradas n * n ...
% Data.Red.Bus.alpha, Data.Red.Bus.cv, Data.Red.Bus.cr, Data.Gen.Pv.qgTop, uTop, Data.Red.Bus.uLow, Data.Red.Bus.pcLow, Data.Red.Bus.qcLow, Data.Red.Bus.qcCap, Data.Gen.Pv.sTop, Data.Gen.Pv.pgTop son vectores columnas de largo n
% Data.Red.Bus.v0, Data.Red.Bus.Q0Top, Data.Red.Bus.Q0Low, Data.Red.Bus.P0Top, Data.Red.Bus.P0Low son escalares
% modSym, findOpt son booleano



	NB = size(Data.Red.Branch.T,1);

	%% Inicializacion

	DistrInfoM.ClRes.v	 = 	replicateMat4DApp(	DistrInfo.ClRes.v	,2);
	DistrInfoM.lambdaT	 = 	replicateMat4DApp(	DistrInfo.lambdaT	,2);
	DistrInfoM.muT	 = 	replicateMat4DApp(	DistrInfo.muT	,2);

	DataM.Red.Bus.uLow	 = 	replicateMat4DApp(	Data.Red.Bus.uLow	,2);
	DataM.Red.Bus.uTop	 = 	replicateMat4DApp(	Data.Red.Bus.uTop	,2);


	pL = Data.Util.pzCnLow(DistrInfo.Bus,1,:,:);
	pT = Data.Util.pzCnTop(DistrInfo.Bus,1,:,:);
	vL = DataM.Red.Bus.uLow(DistrInfo.Bus,1,:,:).^2;
	vT = DataM.Red.Bus.uTop(DistrInfo.Bus,1,:,:).^2;


	%% Modelo programacion matematica

	cvx_begin quiet

	% solo para mosek
        for i = 1: size(Config.SubP,1)
            cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
        end

		%% Declaracion de variables y expresiones
		%Variables generales de la red
		variable pCApp(NB,1,Config.Etapas,2);
		variable qCApp(NB,1,Config.Etapas,2);
		variable poC(NB,1,Config.Etapas,2); % real power demand in i
		variable pzC(NB,1,Config.Etapas,2); % real power demand in i

		% Variables de utilidad
		variable pCn(NB,1,Config.Etapas, 2); % real power demand in i


		% Temperatura Aire Acondicionado
		variable Tvar(NB,1,Config.Etapas);


		%% Declaracion de variables y expresiones


		% Expresion de utilidad
		expression tUtil;

		variable pC(NB,1,Config.Etapas); % real power demand in i
		variable qC(NB,1,Config.Etapas); % reactive power demand in i
		expression TvarAnt(NB,1,Config.Etapas);
		expression pCAC(NB,1,Config.Etapas);
		expression tfopt_expr(Config.Etapas); 
		expression tfopt_virt(Config.Etapas); 
		expression tfopt_conv(Config.Etapas); 
		expression fopt_expr; 

		expression tUtilAC;

		fopt_expr = 0;

		%% Modelado de etapas

		%% Funcion objetivo
		tUtil = -Data.Util.betaT(DistrInfo.Bus,1,:,1).*(pCn(DistrInfo.Bus,1,:,1) - Data.Util.pzCnPref(DistrInfo.Bus,1,:,1)).^2 + Data.Util.aT(DistrInfo.Bus,1,:,1);
		tUtilAC = -Data.St.AC.beta(DistrInfo.Bus,1,:).*(Tvar(DistrInfo.Bus,1,:) - Data.St.AC.tempPref(DistrInfo.Bus,1,:)).^2 + Data.St.AC.a(DistrInfo.Bus,1,:);
		tfopt_expr = - sum(tUtil + tUtilAC);

		tfopt_virt = sum(sum(DistrInfoM.muT(DistrInfo.Bus,1,:,:) .* pCApp(DistrInfo.Bus,1,:,:) + DistrInfoM.lambdaT(DistrInfo.Bus,1,:,:) .* qCApp(DistrInfo.Bus,1,:,:),4),1);
					
		tfopt_conv = 1/(2*DistrInfo.Gama) * ...
			sum(sum(norms(pCApp(DistrInfo.Bus,1,:,:) - DistrInfo.ClRes.pCApp(DistrInfo.Bus,1,:,:),2,2) ...
			+ norms(qCApp(DistrInfo.Bus,1,:,:) - DistrInfo.ClRes.qCApp(DistrInfo.Bus,1,:,:),2,2),4),1);

		fopt_expr = sum(tfopt_expr) + sum(tfopt_virt) + sum(tfopt_conv);

		poC(DistrInfo.Bus,1,:,:) == (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,:,:)).* pCn(DistrInfo.Bus,1,:,:);

		% TODO: si la tension viene dada no hay mac cormick
		pzC(DistrInfo.Bus,1,:,:) >= Data.Red.Bus.alpha(DistrInfo.Bus,1,:,:).*(pL.*DistrInfoM.ClRes.v(DistrInfo.Bus,1,:,:) + pCn(DistrInfo.Bus,1,:,:).*vL - pL.*vL);
		pzC(DistrInfo.Bus,1,:,:) >= Data.Red.Bus.alpha(DistrInfo.Bus,1,:,:).*(pT.*DistrInfoM.ClRes.v(DistrInfo.Bus,1,:,:) + pCn(DistrInfo.Bus,1,:,:).*vT - pT.*vT);
		pzC(DistrInfo.Bus,1,:,:) <= Data.Red.Bus.alpha(DistrInfo.Bus,1,:,:).*(pT.*DistrInfoM.ClRes.v(DistrInfo.Bus,1,:,:) + pCn(DistrInfo.Bus,1,:,:).*vL - pT.*vL);
		pzC(DistrInfo.Bus,1,:,:) <= Data.Red.Bus.alpha(DistrInfo.Bus,1,:,:).*(pL.*DistrInfoM.ClRes.v(DistrInfo.Bus,1,:,:) + pCn(DistrInfo.Bus,1,:,:).*vT - pL.*vT);

		pCn(DistrInfo.Bus,1,:,:) >= pL;
		pCn(DistrInfo.Bus,1,:,:) <= pT;

		pCApp(DistrInfo.Bus,1,:,:) >= pzC(DistrInfo.Bus,1,:,:) + poC(DistrInfo.Bus,1,:,:);
		qCApp(DistrInfo.Bus,1,:,:) == pCApp(DistrInfo.Bus,1,:,:).*Data.Util.tgPhi;


		pC(DistrInfo.Bus,1,:) >= sum(pCApp(DistrInfo.Bus,1,:,:),4);
		qC(DistrInfo.Bus,1,:) >= sum(qCApp(DistrInfo.Bus,1,:,:),4);


		%% Restricciones de Aire Acondicionado
		TvarAnt(:,1,1) = Data.St.AC.tempIni;
		TvarAnt(:,1,(2:Config.Etapas)) = Tvar(:,1,(1:Config.Etapas-1));
		pCAC(:,:,:) = 0;
		pCAC(DistrInfo.Bus,1,:) = pCApp(DistrInfo.Bus,1,:,2);

		Tvar(DistrInfo.Bus,1,:) == TvarAnt(DistrInfo.Bus,1,:) - Data.St.AC.epsilon(DistrInfo.Bus,1,:).*(TvarAnt(DistrInfo.Bus,1,:) - Data.temp(DistrInfo.Bus,1,:))*Data.dt + Data.St.AC.eta(DistrInfo.Bus,1,:).*pCAC(DistrInfo.Bus,1,:)*Data.dt;

		Tvar(DistrInfo.Bus,1,:) >= Data.St.AC.tempLow(DistrInfo.Bus,1,:);
		Tvar(DistrInfo.Bus,1,:) <= Data.St.AC.tempTop(DistrInfo.Bus,1,:);


		if findOpt
			minimize fopt_expr
		end
		
	cvx_end

	Var.ClRes.pC = pC;
	Var.ClRes.qC = qC;
	Var.ClRes.pCApp = pCApp;
	Var.ClRes.qCApp = qCApp;
	Var.ClRes.Tvar = Tvar;

	Var.Red.Bus.pC = pC;
	Var.Red.Bus.qC = qC;
	
	opt(1,1) = sum(tfopt_expr);
	opt(1,2) = sum(tfopt_virt);
	opt(1,3) = sum(tfopt_conv);
	status = cvx_status;

	% for et = 1:Config.Etapas
	
		% for app = 1:2
			% LL = pL(:,1,et,app).*DistrInfo.ClRes.v(DistrInfo.Bus,1,et) + pCn(DistrInfo.Bus,1,et,app).*vL(:,1,et,app) - pL(:,1,et,app).*vL(:,1,et,app);
			% TT = pT(:,1,et,app).*DistrInfo.ClRes.v(DistrInfo.Bus,1,et) + pCn(DistrInfo.Bus,1,et,app).*vT(:,1,et,app) - pT(:,1,et,app).*vT(:,1,et,app);
			% TL = pT(:,1,et,app).*DistrInfo.ClRes.v(DistrInfo.Bus,1,et) + pCn(DistrInfo.Bus,1,et,app).*vL(:,1,et,app) - pT(:,1,et,app).*vL(:,1,et,app);
			% LT = pL(:,1,et,app).*DistrInfo.ClRes.v(DistrInfo.Bus,1,et) + pCn(DistrInfo.Bus,1,et,app).*vT(:,1,et,app) - pL(:,1,et,app).*vT(:,1,et,app);

			% LL = Data.Red.Bus.alpha(DistrInfo.Bus,1,et,app).*(LL) + (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,et,app)).* pCn(DistrInfo.Bus,1,et,app);
			% TT = Data.Red.Bus.alpha(DistrInfo.Bus,1,et,app).*(TT) + (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,et,app)).* pCn(DistrInfo.Bus,1,et,app);
			% TL = Data.Red.Bus.alpha(DistrInfo.Bus,1,et,app).*(TL) + (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,et,app)).* pCn(DistrInfo.Bus,1,et,app);
			% LT = Data.Red.Bus.alpha(DistrInfo.Bus,1,et,app).*(LT) + (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,et,app)).* pCn(DistrInfo.Bus,1,et,app);

			% Var.Red.Bus.errLLR(DistrInfo.Bus,et,app) = abs((pCApp(DistrInfo.Bus,1,et,app) - LL) ./ LL);
			% Var.Red.Bus.errTTR(DistrInfo.Bus,et,app) = abs((pCApp(DistrInfo.Bus,1,et,app) - TT) ./ TT);
			% Var.Red.Bus.errTLR(DistrInfo.Bus,et,app) = abs((pCApp(DistrInfo.Bus,1,et,app) - TL) ./ TL);
			% Var.Red.Bus.errLTR(DistrInfo.Bus,et,app) = abs((pCApp(DistrInfo.Bus,1,et,app) - LT) ./ LT);

			% Var.Red.Bus.errLLA(DistrInfo.Bus,et,app) = abs((pCApp(DistrInfo.Bus,1,et,app) - LL));
			% Var.Red.Bus.errTTA(DistrInfo.Bus,et,app) = abs((pCApp(DistrInfo.Bus,1,et,app) - TT));
			% Var.Red.Bus.errTLA(DistrInfo.Bus,et,app) = abs((pCApp(DistrInfo.Bus,1,et,app) - TL));
			% Var.Red.Bus.errLTA(DistrInfo.Bus,et,app) = abs((pCApp(DistrInfo.Bus,1,et,app) - LT));
		% end
	% end
	
	cvx_clear;
end

function [A] = replicateMat4DApp(V, App)
	A = zeros(size(V,1), size(V,2), size(V,3), App);
	for app = 1:App
		A(:,:,:,app) = V;
	end
end
