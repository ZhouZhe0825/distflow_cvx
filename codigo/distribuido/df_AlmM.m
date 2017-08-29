function [Var, opt, status] = df_AlmM(Data, Config, DistrInfo)
%% Inicializacion
VertI = VertIMat(Data.Red.Branch.T);
VertJ = VertJMat(Data.Red.Branch.T);
OutBr = VertI';
InBr = VertJ';

n = size(VertI,2);

% Baterias
St = DistrInfo.Bus;
nSt = length(St);

M = [];
if nSt > 0
	M = zeros(Config.Etapas,Config.Etapas,nSt);
	for i = 1:nSt
		for et = 1: Config.Etapas-1
			Maux = [Data.St.Bat.m1(St(i),et) -Data.St.Bat.m2(St(i),et)/2; ...
				-Data.St.Bat.m2(St(i),et)/2 Data.St.Bat.m1(St(i),et)];
            M((et:et+1),(et:et+1),i) = M((et:et+1),(et:et+1),i) + Maux;
		end
	end
end

%% Modelo programacion matematica

tic
cvx_begin quiet

	for i = 1: size(Config.SubP,1)
		cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
	end

	cvx_precision high
	%% Declaracion de variables y expresiones
	variable pStb(nSt, Config.Etapas);
	variable qStb(nSt, Config.Etapas);

	expression tfopt_expr(Config.Etapas,1); 
	expression tfopt_mu(Config.Etapas,1); 
	expression tfopt_lambda(Config.Etapas,1); 
	expression tfopt_conv(Config.Etapas,1); 
	expression fopt_expr; 



	%% Restricciones de almacenamiento de bateria
	variable sStb(nSt, Config.Etapas);
	variable pStgb(nSt, Config.Etapas);
	variable pStgbC(nSt, Config.Etapas);
	variable pStgbD(nSt, Config.Etapas);
	variable xiStb(nSt, Config.Etapas);
	variable EStb(n, Config.Etapas);
	variable DlEStb(n, Config.Etapas);

	expression StbNorm(nSt, Config.Etapas,2);
	expression EStbAnt(nSt, Config.Etapas);
	expression cStb(n, Config.Etapas);

	DlEStb <= 0;
	DlEStb <= EStb - Data.St.Bat.ETop.*Data.St.Bat.kapa;


	cStb = Data.St.Bat.wOm + Data.St.Bat.m3.*(DlEStb.^2); % falta termino de m2
	for i = 1:nSt
		cStb(St(i),Config.Etapas) = cStb(St(i),Config.Etapas) + pStgb(i,:) * M(:,:,i) * pStgb(i,:)';
	end

	tfopt_expr = sum(cStb,1) + ...
			sum(Data.St.Bat.beta.* ...
			((Data.St.Bat.EPref - EStb).^2) ...
			+ Data.St.Bat.wU,1);

	tfopt_mu = DistrInfo.lambdaT(St,:) .* qStb;
	tfopt_lambda = DistrInfo.lambdaT(St,:) .* qStb;

	tfopt_conv = 1/(2*DistrInfo.Gama) * ...
			(norms(pStb - DistrInfo.Alm.pStb(St,:),2,2) ...
			+ norms(qStb - DistrInfo.Alm.qStb(St,:),2,2));

	EStbAnt(:,1) = Data.St.Bat.EIni(St,1);
	EStbAnt(:,(2:Config.Etapas)) = EStb(St,(1:Config.Etapas-1));

	pStgb == pStgbC - pStgbD;
	pStgbC >= 0;
	pStgbD >= 0;

	pStb == pStgb - (Data.St.Bat.cv(St,:).*sStb + Data.St.Bat.cr(St,:).*xiStb);
	EStb(St,:) == (1-Data.St.Bat.epsilon(St,:)).*EStbAnt - pStgbD.*Data.St.Bat.etaD(St,:)*Data.dt + Data.St.Bat.etaC(St,:).*pStgbC*Data.dt;

	StbNorm(:,:,1) = pStgb;
	StbNorm(:,:,2) = qStb;
	sStb >= norms(StbNorm,2,3);

	xiStb >= pStgb.^2 + qStb.^2;

	pStb <= Data.St.Bat.pgTop(St,:);
	sStb <= Data.St.Bat.sTop(St,:);
	xiStb <= Data.St.Bat.xiTop(St,:);
	EStb(St,:) <= Data.St.Bat.ETop(St,:);

	pStb >= Data.St.Bat.pgLow(St,:);
	EStb(St,:) >= Data.St.Bat.ELow(St,:);

	fopt_expr = sum(tfopt_expr + tfopt_mu + tfopt_lambda + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT

% Baterias
Var.St.Bat.pStb = zeros(n, Config.Etapas);
Var.St.Bat.pStb(St,:) = pStb;
Var.St.Bat.pStgb = zeros(n, Config.Etapas);
Var.St.Bat.pStgb(St,:) = pStgb;
Var.St.Bat.pStgbC = zeros(n, Config.Etapas);
Var.St.Bat.pStgbC(St,:) = pStgbC;
Var.St.Bat.pStgbD = zeros(n, Config.Etapas);
Var.St.Bat.pStgbD(St,:) = pStgbD;
Var.St.Bat.qStb = zeros(n, Config.Etapas);
Var.St.Bat.qStb(St,:) = qStb;
Var.St.Bat.sStb = zeros(n, Config.Etapas);
Var.St.Bat.sStb(St,:) = sStb;
Var.St.Bat.xiStb = zeros(n, Config.Etapas);
Var.St.Bat.xiStb(St,:) = xiStb;
Var.St.Bat.EStb = zeros(n, Config.Etapas);
Var.St.Bat.EStb(St,:) = EStb(St,:);
Var.St.Bat.cStb = cStb;


Var.Red.Bus.pG = Var.St.Bat.pStb;
Var.Red.Bus.qG = Var.St.Bat.qStb;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_mu);
opt(1,3) = sum(tfopt_lambda);
opt(1,4) = sum(tfopt_conv);


status = cvx_status;
cvx_clear
end