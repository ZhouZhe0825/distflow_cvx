function [Var, opt, status] = df_ClNI(Data, Config, findOpt, DistrInfo)
% function [Var, opt, status] = distflow_nod_eol_util(Data, rhopgSn, rhoP, rhopgWn, rhomqgSn, rhoMqgSn, rhomQ, rhoMQ, rhomqgWn, rhoMqgWn, delta, m, cdv, findOpt, fijarPgSn)
% Data.Red.Branch.T, Data.Red.Branch.r, Data.Red.Branch.x, Data.Gen.Pv.I son matrices cuadradas n * n ...
% Data.Red.Bus.alpha, Data.Red.Bus.cv, Data.Red.Bus.cr, Data.Gen.Pv.qgTop, uTop, Data.Red.Bus.uLow, Data.Red.Bus.pcLow, Data.Red.Bus.qcLow, Data.Red.Bus.qcCap, Data.Gen.Pv.sTop, Data.Gen.Pv.pgTop son vectores columnas de largo n
% Data.Red.Bus.v0, Data.Red.Bus.Q0Top, Data.Red.Bus.Q0Low, Data.Red.Bus.P0Top, Data.Red.Bus.P0Low son escalares
% modSym, findOpt son booleano


    fixed = isfield(Data, 'Fixed');



	NB = size(Data.Red.Branch.T,1);
	n = size(DistrInfo.Bus);

	%% Inicializacion

	pL = Data.Util.pzCnLowE(DistrInfo.Bus,1,:);
	pT = Data.Util.pzCnTopE(DistrInfo.Bus,1,:);
	qL = Data.Util.qzCnLowE(DistrInfo.Bus,1,:);
	qT = Data.Util.qzCnTopE(DistrInfo.Bus,1,:);
	vL = Data.Red.Bus.uLow(DistrInfo.Bus,1,:).^2;
	vT = Data.Red.Bus.uTop(DistrInfo.Bus,1,:).^2;


	%% Modelo programacion matematica

	cvx_begin quiet

	% solo para mosek
        for i = 1: size(Config.SubP,1)
            cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
        end

		%% Declaracion de variables y expresiones
		%Variables generales de la red
		variable pCn(1,1,Config.Etapas); % real power demand in i
		variable poC(1,1,Config.Etapas); % real power demand in i
		variable pzC(1,1,Config.Etapas); % real power demand in i
% 		variable v(Config.Etapas); % module of square complex voltage in i % TODO: preguntar esta variable viene dada por OpSis
		
		% Variables de utilidad
        if fixed
            variable on(1,1,Config.Etapas);
            variable start(1,1,Config.Etapas);
        else
            variable on(1,1,Config.Etapas) binary;
            variable start(1,1,Config.Etapas) binary;
        end
        
            


		% Temperatura Aire Acondicionado
		

		%% Declaracion de variables y expresiones
		%Variables generales de la red
% 		expression tv(n, 1); % module of square complex voltage in i


		% Expresion de utilidad

		expression onNext(1,1,Config.Etapas)

		variable pC(1,1,Config.Etapas); % real power demand in i
		variable qC(1,1,Config.Etapas); % reactive power demand in i
		expression tfopt_expr(Config.Etapas); 
		expression tfopt_virt(Config.Etapas); 
		expression tfopt_conv(Config.Etapas); 
		expression fopt_expr; 


		%% Modelado de etapas
%			 tic
% 			tv = v(:,et);

			%% Funcion objetivo

		tfopt_expr = sum(Data.Util.betaE(DistrInfo.Bus,1,:).*(pCn - Data.Util.pzCnPrefE(DistrInfo.Bus,1,:)).^2 + Data.Util.aE(DistrInfo.Bus,1,:));

		tfopt_virt = sum(DistrInfo.muT(DistrInfo.Bus,1,:) .* pC + DistrInfo.lambdaT(DistrInfo.Bus,1,:) .* qC,1);

		tfopt_conv = 1/(2*DistrInfo.Gama) * ...
			sum(norms(pC - DistrInfo.ClNI.pC(DistrInfo.Bus,1,:),2,2) ...
			+ norms(qC - DistrInfo.ClNI.qC(DistrInfo.Bus,1,:),2,2),1);


		poC == (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,:,1)).* pCn;

		% TODO: si la tension viene dada no hay mac cormick
		pzC >= Data.Red.Bus.alpha(DistrInfo.Bus,1,:,1).*(pL.*DistrInfo.ClNI.v(DistrInfo.Bus,1,:) + pCn.*vL - pL.*vL);
		pzC >= Data.Red.Bus.alpha(DistrInfo.Bus,1,:,1).*(pT.*DistrInfo.ClNI.v(DistrInfo.Bus,1,:) + pCn.*vT - pT.*vT);
		pzC <= Data.Red.Bus.alpha(DistrInfo.Bus,1,:,1).*(pT.*DistrInfo.ClNI.v(DistrInfo.Bus,1,:) + pCn.*vL - pT.*vL);
		pzC <= Data.Red.Bus.alpha(DistrInfo.Bus,1,:,1).*(pL.*DistrInfo.ClNI.v(DistrInfo.Bus,1,:) + pCn.*vT - pL.*vT);
		
		pL.*on <= pCn;
 		pCn <= pT.*on;

		qL.*on <= qC;
 		qC <= qT.*on;

        pC >= pzC + poC;
		qC <= pC.*Data.Util.tgPhi;
		
		sum(start,3) == 1;
		sum(on,3) == Data.ClNI.d(DistrInfo.Bus);
        start(:,1,1) == on(:,1,1);

		onNext(:,1,(1:Config.Etapas-1)) = on(:,1,(2:Config.Etapas));
		onNext(:,1,Config.Etapas) = 0;
		start - onNext + on >= 0;


        fopt_expr = sum(tfopt_expr) + sum(tfopt_virt) + sum(tfopt_conv);
%			 auxiter = toc

        if fixed
            on == Data.Fixed.onCh(DistrInfo.Bus,:,:);
            start == Data.Fixed.stCh(DistrInfo.Bus,:,:); 
        end

		if findOpt
			minimize fopt_expr
		end
		
	cvx_end

	Var.ClNI.pC = zeros(NB,1,Config.Etapas);
	Var.ClNI.qC = zeros(NB,1,Config.Etapas);

	Var.ClNI.pC(DistrInfo.Bus,1,:) = pC;
	Var.ClNI.qC(DistrInfo.Bus,1,:) = qC;

	Var.ClNI.on   = zeros(NB,1,Config.Etapas);
	Var.ClNI.start = zeros(NB,1,Config.Etapas);

	Var.ClNI.on(DistrInfo.Bus,1,:) = on;
	Var.ClNI.start(DistrInfo.Bus,1,:) = start;
	
	Var.Red.Bus.pC = Var.ClNI.pC;
	Var.Red.Bus.qC = Var.ClNI.qC;

	opt(1,1) = sum(tfopt_expr);
	opt(1,2) = sum(tfopt_virt);
	opt(1,3) = sum(tfopt_conv);
	status = cvx_status;
	
	% Var.ClNI.errLLR = zeros(NB,Config.Etapas);
	% Var.ClNI.errTTR = zeros(NB,Config.Etapas);
	% Var.ClNI.errTLR = zeros(NB,Config.Etapas);
	% Var.ClNI.errLTR = zeros(NB,Config.Etapas);

	% Var.ClNI.errLLA = zeros(NB,Config.Etapas);
	% Var.ClNI.errTTA = zeros(NB,Config.Etapas);
	% Var.ClNI.errTLA = zeros(NB,Config.Etapas);
	% Var.ClNI.errLTA = zeros(NB,Config.Etapas);
	
	% for et = 1:Config.Etapas
	
		% LL = pL(:,1,et).*DistrInfo.ClNI.v(DistrInfo.Bus,1,et) + pCn(:,1,et).*vL(:,1,et) - pL(:,1,et).*vL(:,1,et);
		% TT = pT(:,1,et).*DistrInfo.ClNI.v(DistrInfo.Bus,1,et) + pCn(:,1,et).*vT(:,1,et) - pT(:,1,et).*vT(:,1,et);
		% TL = pT(:,1,et).*DistrInfo.ClNI.v(DistrInfo.Bus,1,et) + pCn(:,1,et).*vL(:,1,et) - pT(:,1,et).*vL(:,1,et);
		% LT = pL(:,1,et).*DistrInfo.ClNI.v(DistrInfo.Bus,1,et) + pCn(:,1,et).*vT(:,1,et) - pL(:,1,et).*vT(:,1,et);

		% LL = Data.Red.Bus.alpha(DistrInfo.Bus,1,et,1).*(LL) + (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,et,1)).* pCn(:,1,et);
		% TT = Data.Red.Bus.alpha(DistrInfo.Bus,1,et,1).*(TT) + (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,et,1)).* pCn(:,1,et);
		% TL = Data.Red.Bus.alpha(DistrInfo.Bus,1,et,1).*(TL) + (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,et,1)).* pCn(:,1,et);
		% LT = Data.Red.Bus.alpha(DistrInfo.Bus,1,et,1).*(LT) + (1-Data.Red.Bus.alpha(DistrInfo.Bus,1,et,1)).* pCn(:,1,et);

		% Var.ClNI.errLLR(DistrInfo.Bus,et) = abs((pCn(:,1,et) - LL) ./ LL);
		% Var.ClNI.errTTR(DistrInfo.Bus,et) = abs((pCn(:,1,et) - TT) ./ TT);
		% Var.ClNI.errTLR(DistrInfo.Bus,et) = abs((pCn(:,1,et) - TL) ./ TL);
		% Var.ClNI.errLTR(DistrInfo.Bus,et) = abs((pCn(:,1,et) - LT) ./ LT);

		% Var.ClNI.errLLA(DistrInfo.Bus,et) = abs((pCn(:,1,et) - LL));
		% Var.ClNI.errTTA(DistrInfo.Bus,et) = abs((pCn(:,1,et) - TT));
		% Var.ClNI.errTLA(DistrInfo.Bus,et) = abs((pCn(:,1,et) - TL));
		% Var.ClNI.errLTA(DistrInfo.Bus,et) = abs((pCn(:,1,et) - LT));
	% end

	cvx_clear;
end

