function [Var, opt, status] = df_Tras(Data, Config, findOpt, DistrInfo)
% function [Var, opt, status] = distflow_nod_eol_util(Data, rhopgSn, rhoP, rhopgWn, rhomqgSn, rhoMqgSn, rhomQ, rhoMQ, rhomqgWn, rhoMqgWn, delta, m, cdv, findOpt, fijarPgSn)
% Data.Red.Branch.T, Data.Red.Branch.r, Data.Red.Branch.x, Data.Gen.Pv.I son matrices cuadradas n * n ...
% Data.Red.Bus.alpha, Data.Red.Bus.cv, Data.Red.Bus.cr, Data.Gen.Pv.qgTop, uTop, Data.Red.Bus.uLow, Data.Red.Bus.pcLow, Data.Red.Bus.qcLow, Data.Red.Bus.qcCap, Data.Gen.Pv.sTop, Data.Gen.Pv.pgTop son vectores columnas de largo n
% Data.Red.Bus.v0, Data.Red.Bus.Q0Top, Data.Red.Bus.Q0Low, Data.Red.Bus.P0Top, Data.Red.Bus.P0Low son escalares
% modSym, findOpt son booleano


	n = length(DistrInfo.Bus);
	NB = size(Data.Red.Branch.T,1);

	%% Inicializacion

	%% Modelo programacion matematica

	cvx_begin quiet

	% solo para mosek
        for i = 1: size(Config.SubP,1)
            cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
        end

		%% Declaracion de variables y expresiones
		%Variables generales de la red
		variable PTras(n, 1, Config.Etapas);
		variable QTras(n, 1, Config.Etapas);

		variable cQTras(n, 1, Config.Etapas);

		%% Declaracion de variables y expresiones

		expression tfopt_expr(Config.Etapas); 
		expression tfopt_virt(Config.Etapas); 
		expression tfopt_conv(Config.Etapas); 
		expression fopt_expr; 


		fopt_expr = 0;

		%% Modelado de etapas
		tfopt_conv = 1/(2*DistrInfo.Gama) * ...
				(norms(PTras - sum(DistrInfo.Tras.P(DistrInfo.Bus,:,:),2),2,2) ...
				+ norms(QTras - sum(DistrInfo.Tras.Q(DistrInfo.Bus,:,:),2),2,2));

		tfopt_expr = sum(Data.Cost.piPTras(DistrInfo.Bus,1,:) .* PTras + cQTras);
		tfopt_virt = - DistrInfo.muT(DistrInfo.Bus,1,:) .* PTras(DistrInfo.Bus,1,:) ...
			- DistrInfo.lambdaT(DistrInfo.Bus,1,:) .* QTras(DistrInfo.Bus,1,:);

		fopt_expr = sum(tfopt_expr) + sum(tfopt_virt) + sum(tfopt_conv);

		cQTras >= - Data.Cost.piQmtras(DistrInfo.Bus,1,:) .* QTras;
		cQTras >= Data.Cost.piQMtras(DistrInfo.Bus,1,:) .* QTras;

		% Restricciones de nodo 0
		PTras >= Data.Red.Bus.P0Low(DistrInfo.Bus,:,:);
		PTras <= Data.Red.Bus.P0Top(DistrInfo.Bus,:,:);
		QTras >= Data.Red.Bus.Q0Low(DistrInfo.Bus,:,:);
		QTras <= Data.Red.Bus.Q0Top(DistrInfo.Bus,:,:);

		if findOpt
			minimize fopt_expr
		end
		
	cvx_end

	Var.Red.Bus.PTras = zeros(NB,1,Config.Etapas);
	Var.Red.Bus.QTras = zeros(NB,1,Config.Etapas);

	Var.Red.Bus.PTras(DistrInfo.Bus,1,:) = PTras;
	Var.Red.Bus.QTras(DistrInfo.Bus,1,:) = QTras;
	Var.Cost.optEt = tfopt_expr;

	Var.Red.Bus.pG = Var.Red.Bus.PTras;
	Var.Red.Bus.qG = Var.Red.Bus.QTras;
	
	opt(1,1) = sum(tfopt_expr);
	opt(1,2) = sum(tfopt_virt);
	opt(1,3) = sum(tfopt_conv);
	status = cvx_status;

	
	cvx_clear;
end
