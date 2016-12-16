function [Var, opt, status] = df_Alm(Data, Config, findOpt, DistrInfo)
% function [Var, opt, status] = distflow_nod_eol_util(Data, rhopgSn, rhoP, rhopgWn, rhomqgSn, rhoMqgSn, rhomQ, rhoMQ, rhomqgWn, rhoMqgWn, delta, m, cdv, findOpt, fijarPgSn)
% Data.Red.Branch.T, Data.Red.Branch.r, Data.Red.Branch.x, Data.Gen.Pv.I son matrices cuadradas n * n ...
% Data.Red.Bus.alpha, Data.Red.Bus.cv, Data.Red.Bus.cr, Data.Gen.Pv.qgTop, uTop, Data.Red.Bus.uLow, Data.Red.Bus.pcLow, Data.Red.Bus.qcLow, Data.Red.Bus.qcCap, Data.Gen.Pv.sTop, Data.Gen.Pv.pgTop son vectores columnas de largo n
% Data.Red.Bus.v0, Data.Red.Bus.Q0Top, Data.Red.Bus.Q0Low, Data.Red.Bus.P0Top, Data.Red.Bus.P0Low son escalares
% modSym, findOpt son booleano



	%% Inicializacion
	NB = size(Data.Red.Branch.T,1);
	M = zeros(2,2,Config.Etapas);
	for et = 1: Config.Etapas
		M(:,:,et) = [Data.St.Bat.m1(DistrInfo.Bus,1,et) -Data.St.Bat.m2(DistrInfo.Bus,1,et)/2; ...
			-Data.St.Bat.m2(DistrInfo.Bus,1,et)/2 Data.St.Bat.m1(DistrInfo.Bus,1,et)];
	end

	%% Modelo programacion matematica

	cvx_begin quiet

	% solo para mosek
        for i = 1: size(Config.SubP,1)
            cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
        end

		%% Declaracion de variables y expresiones
		% Almacenamiento en baterias
		variables pStgb(1,1,Config.Etapas);
		variables pStb(1,1,Config.Etapas);
		variables qStb(1,1,Config.Etapas);
		variables sStb(1,1,Config.Etapas);
		variables xiStb(1,1,Config.Etapas);
		variables EStb(1,1,Config.Etapas);
		variables DlEStb(1,1,Config.Etapas);

		expression tfopt_expr(Config.Etapas); 
		expression tfopt_virt(Config.Etapas); 
		expression tfopt_conv(Config.Etapas); 
		expression fopt_expr; 


		expression cStom(Config.Etapas);
		expression cStp(Config.Etapas);
		
		expression EStbAnt(1,1,Config.Etapas);

		%% Modelado de etapas
		cStom = Data.St.Bat.wOm(DistrInfo.Bus,1,:) + Data.St.Bat.m3(DistrInfo.Bus,1,:).*(DlEStb.^2); % falta termino de m2

		for et = 1: Config.Etapas
			if et == 1
				cStom(et) = cStom(et) + [pStgb(1,1,et) 0] * M(:,:,et) * [pStgb(1,1,et); 0];
			else
				cStom(et) = cStom(et) + [pStgb(1,1,et) pStgb(1,1,et-1)] * M(:,:,et) * [pStgb(1,1,et); pStgb(1,1,et-1)];
			end
		end

		cStp = -Data.St.Bat.beta(DistrInfo.Bus,1,:).*((Data.St.Bat.ETop(DistrInfo.Bus,1,:) - EStb(1,1,Config.Etapas)*Data.St.Bat.gama).^2) + Data.St.Bat.wU(DistrInfo.Bus,1,:);

		tfopt_expr = sum(cStom) - sum(cStp);
		tfopt_virt = DistrInfo.muT(DistrInfo.Bus,1,:) .* pStb + DistrInfo.lambdaT(DistrInfo.Bus,1,:) .* qStb;

		tfopt_conv = 1/(2*DistrInfo.Gama) * ...
				(norms(pStb - DistrInfo.Alm.pStb(DistrInfo.Bus,1,:),2,2) ...
				+ norms(qStb - DistrInfo.Alm.qStb(DistrInfo.Bus,1,:),2,2));

		fopt_expr = sum(tfopt_expr) + sum(tfopt_virt) + sum(tfopt_conv);

		pStb == pStgb - (Data.St.Bat.cv(DistrInfo.Bus,1,:).*sStb + Data.St.Bat.cr(DistrInfo.Bus,1,:).*xiStb);

		EStbAnt(1,1,1) = Data.St.Bat.EIni(DistrInfo.Bus,1,1);
		EStbAnt(1,1,(2:Config.Etapas)) = EStb(1,1,(1:Config.Etapas-1));

		EStb == (1-Data.St.Bat.epsilon(DistrInfo.Bus,1,:)).*EStbAnt - Data.St.Bat.eta(DistrInfo.Bus,1,:).*pStgb*Data.dt;


		DlEStb <= 0;
		DlEStb <= EStb - Data.St.Bat.ETop(DistrInfo.Bus,1,:)*Data.St.Bat.kapa;
			
			%% Restricciones de storage bateria
		sStb >= norms([pStgb qStb], 2, 2);
		xiStb >= pStgb.^2 + qStb.^2;

		pStb <= Data.St.Bat.pgTop(DistrInfo.Bus,1,:);
		sStb <= Data.St.Bat.sTop(DistrInfo.Bus,1,:);
		EStb <= Data.St.Bat.ETop(DistrInfo.Bus,1,:);

		pStb >= Data.St.Bat.pgLow(DistrInfo.Bus,1,:);
		EStb >= Data.St.Bat.ELow(DistrInfo.Bus,1,:);

		if findOpt
			minimize fopt_expr
		end
		
	cvx_end

	Var.St.Bat.pStb = zeros(NB, 1, Config.Etapas);
	Var.St.Bat.pStgb = zeros(NB, 1, Config.Etapas);
	Var.St.Bat.qStb = zeros(NB, 1, Config.Etapas);
	Var.St.Bat.sStb = zeros(NB, 1, Config.Etapas);
	Var.St.Bat.xiStb = zeros(NB, 1, Config.Etapas);
	Var.St.Bat.EStb = zeros(NB, 1, Config.Etapas);

	Var.St.Bat.pStb(DistrInfo.Bus,1,:) = pStb;
	Var.St.Bat.pStgb(DistrInfo.Bus,1,:) = pStgb;
	Var.St.Bat.qStb(DistrInfo.Bus,1,:) = qStb;
	Var.St.Bat.sStb(DistrInfo.Bus,1,:) = sStb;
	Var.St.Bat.xiStb(DistrInfo.Bus,1,:) = xiStb;
	Var.St.Bat.EStb(DistrInfo.Bus,1,:) = EStb;

	Var.Red.Bus.pG = Var.St.Bat.pStb;
	Var.Red.Bus.qG = Var.St.Bat.qStb;

	opt(1,1) = sum(tfopt_expr);
	opt(1,2) = sum(tfopt_virt);
	opt(1,3) = sum(tfopt_conv);
	status = cvx_status;

	cvx_clear;
end
