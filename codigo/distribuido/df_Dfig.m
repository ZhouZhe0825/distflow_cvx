function [Var, opt, status] = df_Dfig(Data, Config, findOpt, DistrInfo)
% function [Var, opt, status] = distflow_nod_eol_util(Data, rhopgSn, rhoP, rhopgWn, rhomqgSn, rhoMqgSn, rhomQ, rhoMQ, rhomqgWn, rhoMqgWn, delta, m, cdv, findOpt, fijarPgSn)
% Data.Red.Branch.T, Data.Red.Branch.r, Data.Red.Branch.x, Data.Gen.Pv.I son matrices cuadradas n * n ...
% Data.Red.Bus.alpha, Data.Red.Bus.cv, Data.Red.Bus.cr, Data.Gen.Pv.qgTop, uTop, Data.Red.Bus.uLow, Data.Red.Bus.pcLow, Data.Red.Bus.qcLow, Data.Red.Bus.qcCap, Data.Gen.Pv.sTop, Data.Gen.Pv.pgTop son vectores columnas de largo n
% Data.Red.Bus.v0, Data.Red.Bus.Q0Top, Data.Red.Bus.Q0Low, Data.Red.Bus.P0Top, Data.Red.Bus.P0Low son escalares
% modSym, findOpt son booleano



%% Inicializacion
	NB = size(Data.Red.Branch.T,1);

	% Inicializacion Eolico
	isSymT = isequal(matOverTime(Data.Red.Branch.T), matOverTime(Data.Red.Branch.T)');
	wnd = DistrInfo.Dfig.ind(DistrInfo.Bus);
	lenWnd = length(find(matOverTime(Data.Gen.DFIG.I) == 1));

	if isSymT
		%% Modelo programacion matematica
		cvx_begin quiet

		% solo para mosek
            for i = 1: size(Config.SubP,1)
                cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
            end

			%% Declaracion de variables y expresiones

			% Variables de generadores eolico
			variable PDFIG(5,5, Config.Etapas);
			variable lDFIG(5,5, Config.Etapas);
			variable pWigDFIG(5,1, Config.Etapas);
			variable QDFIG(5,5, Config.Etapas);
			variable qWigDFIG(5,1, Config.Etapas);
			variable pCDFIG(5,1, Config.Etapas);
			variable qCDFIG(5,1, Config.Etapas);
			variable vDFIG(5,1, Config.Etapas);
			variable sDFIG(5,1, Config.Etapas);
			variable xiDFIG(5,1, Config.Etapas);

			variable cqWi(1,1,Config.Etapas);

			% Expresion para p y q generados que entran en la red general
			variable pWi(1,1,Config.Etapas);
			variable qWi(1,1,Config.Etapas);


			% i -> e  |-> r
			%  |-> f  o

			% 1 -> 2  |-> 5
			%  |-> 3  4

			expression lQoL(5,5,Config.Etapas,3);

			expression tfopt_expr(Config.Etapas); 
			expression tfopt_virt(Config.Etapas); 
			expression tfopt_conv(Config.Etapas); 
			expression fopt_expr; 
			fopt_expr = 0;



			%% Funcion objetivo
			tfopt_expr = sum(Data.Cost.rhopWi(DistrInfo.Bus,1,:) .* pWi) + sum(cqWi);

			cqWi >= - Data.Cost.rhomqWi(DistrInfo.Bus,1,:) .* qWi;
			cqWi >= Data.Cost.rhoMqWi(DistrInfo.Bus,1,:) .* qWi;

			tfopt_conv = 1/(2*DistrInfo.Gama) * ...
					sum(norms(pWi(1,1,:) - DistrInfo.Dfig.pWi(DistrInfo.Bus,1,:),2,2) ...
					+ norms(qWi(1,1,:) - DistrInfo.Dfig.qWi(DistrInfo.Bus,1,:),2,2),1);

			tfopt_virt = sum(DistrInfo.muT(DistrInfo.Bus,1,:) .* pWi(1,1,:) + DistrInfo.lambdaT(DistrInfo.Bus,1,:) .* qWi(1,1,:),1);

			fopt_expr = sum(tfopt_expr) + sum(tfopt_virt) + sum(tfopt_conv);

			%% Restricciones de Red
			pWi == - PDFIG(1,2,:) - PDFIG(1,3,:);
			qWi == - QDFIG(1,2,:) - QDFIG(1,3,:);

			%% Restricciones de generadores Eolico
			% Modelo de Red interna
			vDFIG(1,1,:) == DistrInfo.Dfig.v(DistrInfo.Bus,1,:); % TODO: poner bien son datos que vienen por otros problemas

			PDFIG(1,2,:) == Data.Gen.DFIG.r(1,2,:,wnd) .* lDFIG(1,2,:) - pWigDFIG(2,1,:);
			PDFIG(1,3,:) == Data.Gen.DFIG.r(1,3,:,wnd) .* lDFIG(1,3,:) + pCDFIG(3,1,:);
			PDFIG(4,5,:) == Data.Gen.DFIG.r(4,5,:,wnd) .* lDFIG(4,5,:) - pWigDFIG(5,1,:);

			QDFIG(1,2,:) == Data.Gen.DFIG.x(1,2,:,wnd) .* lDFIG(1,2,:) - qWigDFIG(2,1,:);
			QDFIG(1,3,:) == Data.Gen.DFIG.x(1,3,:,wnd) .* lDFIG(1,3,:) - qCDFIG(3,1,:);
			QDFIG(4,5,:) == Data.Gen.DFIG.x(4,5,:,wnd) .* (Data.Gen.DFIG.n_(DistrInfo.Bus,1,:) .* lDFIG(4,5,:)) - qWigDFIG(5,1,:);

			vDFIG(2,1,:) == vDFIG(1,1,:) - 2*(Data.Gen.DFIG.r(1,2,:,wnd) .* PDFIG(1,2,:) + Data.Gen.DFIG.x(1,2,:,wnd) .* QDFIG(1,2,:)) + (Data.Gen.DFIG.r(1,2,:,wnd).^2 + Data.Gen.DFIG.x(1,2,:,wnd).^2) .* lDFIG(1,2,:);
			vDFIG(3,1,:) == vDFIG(1,1,:) - 2*(Data.Gen.DFIG.r(1,3,:,wnd) .* PDFIG(1,3,:) + Data.Gen.DFIG.x(1,3,:,wnd) .* QDFIG(1,3,:)) + (Data.Gen.DFIG.r(1,3,:,wnd).^2 + Data.Gen.DFIG.x(1,3,:,wnd).^2) .* lDFIG(1,3,:);
			vDFIG(5,1,:) == vDFIG(4,1,:) - 2*(Data.Gen.DFIG.r(4,5,:,wnd) .* PDFIG(4,5,:) + Data.Gen.DFIG.x(4,5,:,wnd) .* (Data.Gen.DFIG.n_(DistrInfo.Bus,1,:)  .* QDFIG(4,5,:))) + (Data.Gen.DFIG.r(4,5,:,wnd).^2 + Data.Gen.DFIG.n_(DistrInfo.Bus,1,:).^2 .* Data.Gen.DFIG.x(4,5,:,wnd).^2) .* lDFIG(4,5,:);

			% Corriente
			lQoL(:,:,:,1) = 2*PDFIG.*Data.Gen.DFIG.Tg;
			lQoL(:,:,:,2) = 2*QDFIG.*Data.Gen.DFIG.Tg;
			lQol(:,:,:,3) = (lDFIG - repmat(vDFIG(:,1,:), [1 5 1])).*Data.Gen.DFIG.Tg;
			norms(lQol,2,4) - (lDFIG + repmat(vDFIG(:,1,:), [1 5 1])).*Data.Gen.DFIG.Tg <= 0;

			Data.Gen.DFIG.lTop(1,3,:,wnd) >= lDFIG(1,3,:);
			Data.Gen.DFIG.lTop(4,5,:,wnd) >= lDFIG(4,5,:);

			(Data.Gen.DFIG.lTop(1,2,:,wnd) - lDFIG(1,2,:)).*abs(sign(Data.Gen.DFIG.P_mec(DistrInfo.Bus,1,:))) >= 0;
			(Data.Gen.DFIG.sTop(3,1,:,wnd) - sDFIG(3,1,:)).*abs(sign(Data.Gen.DFIG.P_mec(DistrInfo.Bus,1,:))) >= 0;
			(Data.Gen.DFIG.sTop(5,1,:,wnd) - sDFIG(5,1,:)).*abs(sign(Data.Gen.DFIG.P_mec(DistrInfo.Bus,1,:))) >= 0;

			vDFIG(2,1,:) >= Data.Gen.DFIG.uLow(2,1,:,wnd).^2;
			vDFIG(2,1,:) <= Data.Gen.DFIG.uTop(2,1,:,wnd).^2;

			vDFIG(3,1,:) >= Data.Gen.DFIG.uLow(3,1,:,wnd).^2;
			vDFIG(3,1,:) <= Data.Gen.DFIG.uTop(3,1,:,wnd).^2;

			Data.Gen.DFIG.PQnorm(1,2,:,wnd) >= norms([PDFIG(1,2,:) QDFIG(1,2,:)],2 ,2);
			Data.Gen.DFIG.PQnorm(1,3,:,wnd) >= norms([PDFIG(1,3,:) QDFIG(1,3,:)],2 ,2);

			sDFIG(3,1,:) >= norms([pCDFIG(3,1,:) qCDFIG(3,1,:)],2 ,2);
			xiDFIG(3,1,:) >= pCDFIG(3,1,:).^2 + qCDFIG(3,1,:).^2;

			sDFIG(5,1,:) >= norms([PDFIG(4,5,:) QDFIG(4,5,:)],2 ,2);
			xiDFIG(5,1,:) >= PDFIG(4,5,:).^2 + QDFIG(4,5,:).^2;

			pCDFIG(3,1,:) == PDFIG(4,5,:) ...
				+ (Data.Gen.DFIG.cv(5,1,:,wnd) .* sDFIG(5,1,:) + Data.Gen.DFIG.cr(5,1,:,wnd) .* xiDFIG(5,1,:)) ...
				+ (Data.Gen.DFIG.cv(3,1,:,wnd) .* sDFIG(3,1,:) + Data.Gen.DFIG.cr(3,1,:,wnd) .* xiDFIG(3,1,:));
			vDFIG(5,1,:) == Data.Gen.DFIG.N_er^2 * (Data.Gen.DFIG.n_(DistrInfo.Bus,1,:).^2) .* vDFIG(2,1,:);
			pWigDFIG(2,1,:) == Data.Gen.DFIG.P_mec(DistrInfo.Bus,1,:) ./ (1-Data.Gen.DFIG.n_(DistrInfo.Bus,1,:)); %TODO es igual
			pWigDFIG(5,1,:) == -Data.Gen.DFIG.n_(DistrInfo.Bus,1,:) .* pWigDFIG(2,1,:);
			qWigDFIG(5,1,:) == Data.Gen.DFIG.n_(DistrInfo.Bus,1,:) .* (qWigDFIG(2,1,:));

			(1-Data.Gen.DFIG.Tg).*PDFIG == 0;
			(1-Data.Gen.DFIG.Tg).*QDFIG == 0;
			(1-Data.Gen.DFIG.Tg).*lDFIG == 0;
			
			if findOpt
				minimize fopt_expr
			end

		cvx_end

		Var.Gen.Dfig.pWi = zeros(NB,1,Config.Etapas);
		Var.Gen.Dfig.qWi = zeros(NB,1,Config.Etapas);

		Var.Gen.Dfig.pWi(DistrInfo.Bus,1,:) = pWi;
		Var.Gen.Dfig.qWi(DistrInfo.Bus,1,:) = qWi;

		Var.Gen.Dfig.Branch.P = zeros(5,5,Config.Etapas,lenWnd);
		Var.Gen.Dfig.Branch.Q = zeros(5,5,Config.Etapas,lenWnd);
		Var.Gen.Dfig.Branch.l = zeros(5,5,Config.Etapas,lenWnd);
		Var.Gen.Dfig.Bus.v = zeros(5,1,Config.Etapas,lenWnd);
		Var.Gen.Dfig.Bus.pC = zeros(5,1,Config.Etapas,lenWnd);
		Var.Gen.Dfig.Bus.qC = zeros(5,1,Config.Etapas,lenWnd);
		Var.Gen.Dfig.Bus.pg = zeros(5,1,Config.Etapas,lenWnd);
		Var.Gen.Dfig.Bus.qg = zeros(5,1,Config.Etapas,lenWnd);
		Var.Gen.Dfig.Bus.s = zeros(5,1,Config.Etapas,lenWnd);
		Var.Gen.Dfig.Bus.xi = zeros(5,1,Config.Etapas,lenWnd);

		Var.Gen.Dfig.Branch.P(:,:,:,wnd) = PDFIG;
		Var.Gen.Dfig.Branch.Q(:,:,:,wnd) = QDFIG;
		Var.Gen.Dfig.Branch.l(:,:,:,wnd) = lDFIG;
		Var.Gen.Dfig.Bus.v(:,:,:,wnd) = vDFIG;
		Var.Gen.Dfig.Bus.pC(:,:,:,wnd) = pCDFIG;
		Var.Gen.Dfig.Bus.qC(:,:,:,wnd) = qCDFIG;
		Var.Gen.Dfig.Bus.pg(:,:,:,wnd) = pWigDFIG;
		Var.Gen.Dfig.Bus.qg(:,:,:,wnd) = qWigDFIG;
		Var.Gen.Dfig.Bus.s(:,:,:,wnd) = sDFIG;
		Var.Gen.Dfig.Bus.xi(:,:,:,wnd) = xiDFIG;

		Var.Red.Bus.pG = Var.Gen.Dfig.pWi;
		Var.Red.Bus.qG = Var.Gen.Dfig.qWi;

		opt(1,1) = sum(tfopt_expr);
		opt(1,2) = sum(tfopt_virt);
		opt(1,3) = sum(tfopt_conv);
		status = cvx_status;

		% Var.Gen.Dfig.Branch.errLR = Var.Gen.Dfig.Branch.l * 0;
		% Var.Gen.Dfig.Branch.errLA = Var.Gen.Dfig.Branch.errLR;

		% Var.Gen.Dfig.Bus.errSR = Var.Gen.Dfig.Bus.s * 0;
		% Var.Gen.Dfig.Bus.errSA = Var.Gen.Dfig.Bus.errSR;
		% Var.Gen.Dfig.Bus.errXiR = Var.Gen.Dfig.Bus.errSR;
		% Var.Gen.Dfig.Bus.errXiA = Var.Gen.Dfig.Bus.errXiR;

		% auxS = Var.Gen.Dfig.Bus.errSR;
		% auxS(3,1,:) = Var.Gen.Dfig.Bus.pC(3,1,:).^2 + Var.Gen.Dfig.Bus.qC(3,1,:).^2;
		% auxS(5,1,:) = Var.Gen.Dfig.Branch.P(4,5,:).^2 + Var.Gen.Dfig.Branch.Q(4,5,:).^2;
		% Var.Gen.Dfig.Bus.errSR = abs(Var.Gen.Dfig.Bus.s - sqrt(auxS)) ./ sqrt(auxS);

		% Var.Gen.Dfig.Bus.errXiR = abs(Var.Gen.Dfig.Bus.xi - auxS) ./ auxS;

% 		for i = 1:Config.Etapas
% 			PQV = (Var.Gen.Dfig.Branch.P(:,:,i).^2 + Var.Gen.Dfig.Branch.Q(:,:,i).^2) ./ (Var.Gen.Dfig.Bus.v(:,:,i)*ones(1,length(Data.Gen.DFIG.Tg(:,:,i))).*Data.Gen.DFIG.Tg(:,:,i));
% 			Var.Gen.Dfig.Branch.errLR(:,:,i) = abs(Var.Gen.Dfig.Branch.l(:,:,i) - PQV)./PQV;
% 			Var.Gen.Dfig.Branch.errLA(:,:,i) = Data.Gen.DFIG.r .* abs(Var.Gen.Dfig.Branch.l(:,:,i) - PQV);
% 			Var.Gen.Dfig.Bus.errSA(:,:,i) = Data.Gen.DFIG.cv .* abs(Var.Gen.Dfig.Bus.s(:,:,i) - sqrt(auxS(:,:,i)));
% 			Var.Gen.Dfig.Bus.errXiA(:,:,i) = Data.Gen.DFIG.cr .* abs(Var.Gen.Dfig.Bus.xi(:,:,i) - auxS(:,:,i));
% 		end


		cvx_clear;
	end;
end
