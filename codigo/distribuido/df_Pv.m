function [Var, opt, status] = df_Pv(Data, Config, findOpt, DistrInfo)
% function [Var, opt, status] = distflow_nod_eol_util(Data, rhopgSn, rhoP, rhopgWn, rhomqgSn, rhoMqgSn, rhomQ, rhoMQ, rhomqgWn, rhoMqgWn, delta, m, cdv, findOpt, fijarPgSn)
% Data.Red.Branch.T, Data.Red.Branch.r, Data.Red.Branch.x, Data.Gen.Pv.I son matrices cuadradas n * n ...
% Data.Red.Bus.alpha, Data.Red.Bus.cv, Data.Red.Bus.cr, Data.Gen.Pv.qgTop, uTop, Data.Red.Bus.uLow, Data.Red.Bus.pcLow, Data.Red.Bus.qcLow, Data.Red.Bus.qcCap, Data.Gen.Pv.sTop, Data.Gen.Pv.pgTop son vectores columnas de largo n
% Data.Red.Bus.v0, Data.Red.Bus.Q0Top, Data.Red.Bus.Q0Low, Data.Red.Bus.P0Top, Data.Red.Bus.P0Low son escalares
% modSym, findOpt son booleano




	NB = size(Data.Red.Branch.T,1);

%% Inicializacion

	%% Modelo programacion matematica

	cvx_begin quiet

	% solo para mosek
        for i = 1: size(Config.SubP,1)
            cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
        end

		%% Declaracion de variables y expresiones

		% Variables de generadores solares
		variable pPv(1,1,Config.Etapas); % real power generation in i
		variable qPv(1,1,Config.Etapas); % reactive power generation in i
		variable sPv(1,1,Config.Etapas);
		variable xiPv(1,1,Config.Etapas); % module of square complex current in i
		variable cqPv(1,1,Config.Etapas);

		% Variables de generadores solares
		expression tfopt_expr(Config.Etapas); 
		expression tfopt_virt(Config.Etapas); 
		expression tfopt_conv(Config.Etapas); 
		expression fopt_expr; 

		fopt_expr = 0;

		%% Modelado de etapas
		tfopt_expr = sum(Data.Cost.rhopPv(DistrInfo.Bus,1,:) .* pPv) + sum(cqPv);
		tfopt_conv = 1/(2*DistrInfo.Gama) * ...
				(norms(pPv - DistrInfo.PV.pPv(DistrInfo.Bus,1,:),2,2) ...
				+ norms(qPv - DistrInfo.PV.qPv(DistrInfo.Bus,1,:),2,2));
%				 + 1/2*sum(norm(DistrInfo.mu(DistrInfo.Bus) - DistrInfo.muAnt(DistrInfo.Bus),2)) ...
%				 + 1/2*sum(norm(DistrInfo.lambda(DistrInfo.Bus) - DistrInfo.lambdaAnt(DistrInfo.Bus),2));
		tfopt_virt = DistrInfo.muT(DistrInfo.Bus,1,:) .* pPv + DistrInfo.lambdaT(DistrInfo.Bus,1,:) .* qPv;
				


		cqPv >= - Data.Cost.rhomqPv(DistrInfo.Bus,1,:) .* qPv;
		cqPv >= Data.Cost.rhoMqPv(DistrInfo.Bus,1,:) .* qPv;

			%% Restricciones de generadores Solares
		pPv == Data.Gen.Pv.pPvg(DistrInfo.Bus,1,:) - (Data.Gen.Pv.cv(DistrInfo.Bus,1,:).*sPv + Data.Gen.Pv.cr(DistrInfo.Bus,1,:).*xiPv);
		sPv >= norms([pPv qPv], 2, 2);
		xiPv >= pPv.^2 + qPv.^2;
		sPv <= Data.Gen.Pv.sTop(DistrInfo.Bus,1,:).*abs(sign(Data.Gen.Pv.pPvg(DistrInfo.Bus,1,:)));


		fopt_expr = tfopt_expr + sum(tfopt_virt) + sum(tfopt_conv);

		if findOpt
			minimize fopt_expr
		end
		

	cvx_end


% 	Var.Gen.Pv.pPv = zeros(NB,Config.Etapas);
% 	Var.Gen.Pv.qPv = zeros(NB,Config.Etapas);
% 	Var.Gen.Pv.s = zeros(NB,Config.Etapas);
% 	Var.Gen.Pv.xi = zeros(NB,Config.Etapas);
%	 
% 	Var.Gen.Pv.pPv(DistrInfo.Bus,:) = pPv;
% 	Var.Gen.Pv.qPv(DistrInfo.Bus,:) = qPv;
% 	Var.Gen.Pv.s(DistrInfo.Bus,:) = sPv;
% 	Var.Gen.Pv.xi(DistrInfo.Bus,:) = xiPv;
% 
% 	Var.Cost.CpgS = Data.Cost.rhopPv' * Var.Gen.Pv.pPv;
% 	Var.Cost.CqgS = sum(cqPv);
% 	Var.Cost.optEt = tfopt_expr;

	
	Var.Gen.Pv.pPv = zeros(NB,1,Config.Etapas);
	Var.Gen.Pv.qPv = zeros(NB,1,Config.Etapas);
	Var.Gen.Pv.s = zeros(NB,1,Config.Etapas);
	Var.Gen.Pv.xi = zeros(NB,1,Config.Etapas);

	Var.Gen.Pv.pPv(DistrInfo.Bus,1,:) = pPv;
	Var.Gen.Pv.qPv(DistrInfo.Bus,1,:) = qPv;
	Var.Gen.Pv.s(DistrInfo.Bus,1,:) = sPv;
	Var.Gen.Pv.xi(DistrInfo.Bus,1,:) = xiPv;

% 	Var.Gen.Pv.CpgS = Data.Cost.rhopPv' * Var.Gen.Pv.pPv; %TODO arreglar para que ande
	Var.Gen.Pv.CqgS = sum(cqPv);
	% Var.Gen.Pv.optEt = tfopt_expr;

	Var.Red.Bus.pG = Var.Gen.Pv.pPv;
	Var.Red.Bus.qG = Var.Gen.Pv.qPv;

	opt(1,1) = sum(tfopt_expr);
	opt(1,2) = sum(tfopt_virt);
	opt(1,3) = sum(tfopt_conv);
	status = cvx_status;

	% auxS = Var.Gen.Pv.pPv.^2 + Var.Gen.Pv.qPv.^2;
	% Var.Gen.Pv.errSR = abs(Var.Gen.Pv.s - sqrt(auxS)) ./ sqrt(auxS);
%	 Var.Gen.Pv.errSA = Data.Gen.Pv.cv * abs(Var.Gen.Pv.s - sqrt(auxS)); %TODO arreglar para que ande

	% Var.Gen.Pv.errXi = abs(Var.Gen.Pv.s - auxS) ./ auxS;
%	 Var.Gen.Pv.errXi = Data.Gen.Pv.cr * abs(Var.Gen.Pv.s - auxS); %TODO arreglar para que ande

	cvx_clear;
end
