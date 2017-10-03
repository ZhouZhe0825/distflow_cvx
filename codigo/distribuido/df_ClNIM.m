function [Var, opt, status] = df_ClNIM(Data, Config, DistrInfo)

%% Inicializacion
VertI = VertIMat(Data.Red.Branch.T);
VertJ = VertJMat(Data.Red.Branch.T);
OutBr = VertI';
InBr = VertJ';

n = size(VertI,2);
m = size(VertI,1);

% Cargas No interrumpibles
ClNI = DistrInfo.Bus;
nClNI = length(ClNI);

pLClNI = Data.Util.pzCnLowE(ClNI,:);
pTClNI = Data.Util.pzCnTopE(ClNI,:);
qLClNI = Data.Util.qzCnLowE(ClNI,:);
qTClNI = Data.Util.qzCnTopE(ClNI,:);
vLClNI = Data.Red.Bus.uLow(ClNI,:).^2;
vTClNI = Data.Red.Bus.uTop(ClNI,:).^2;

%% Modelo programacion matematica
tic
cvx_begin quiet

	for i = 1: size(Config.SubP,1)
		cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
	end

	cvx_precision high
	%% Declaracion de variables y expresiones
	%Variables generales de la red

	variable pCClNI(nClNI, Config.Etapas);
	variable qCClNI(nClNI, Config.Etapas);

	expression tfopt_expr(Config.Etapas,1); 
	expression tfopt_mu(Config.Etapas,1); 
	expression tfopt_lambda(Config.Etapas,1); 
	expression tfopt_conv(Config.Etapas,1); 
	expression fopt_expr; 

	%% Restricciones de cargas no interrumpibles
	variable onClNI(nClNI,Config.Etapas) binary;
	variable stClNI(nClNI,Config.Etapas) binary;

	expression onNext(nClNI,Config.Etapas)

	pCClNI == Data.ClNI.pC(ClNI,:).*onClNI;
	qCClNI == Data.ClNI.qC(ClNI,:).*onClNI;

	onNext(:,(1:Config.Etapas-1)) = onClNI(:,(2:Config.Etapas));
	onNext(:,Config.Etapas) = 0;
	stClNI - onNext + onClNI >= 0;

	sum(stClNI,2) == 1;
	sum(onClNI,2) == Data.ClNI.d(ClNI);
	stClNI(:,1) == onClNI(:,1);

	tfopt_expr = sum(Data.Util.betaE(ClNI,:).*((pCClNI - Data.Util.pzCnPrefE(ClNI,:)).^2),1);
	tfopt_mu = sum(DistrInfo.muT(ClNI,:) .* pCClNI,1);
	tfopt_lambda = sum(DistrInfo.lambdaT(ClNI,:) .* qCClNI,1);

	tfopt_conv = ...
		1/(2*DistrInfo.Gama_m) * norms(pCClNI - DistrInfo.ClNI.pC(ClNI,:),2,1) ...
		+ 1/(2*DistrInfo.Gama_l) * norms(qCClNI - DistrInfo.ClNI.qC(ClNI,:),2,1);



	fopt_expr = sum(tfopt_expr + tfopt_mu + tfopt_lambda + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT
% Cargas No interrumpibles
if nClNI > 0
	Var.ClNI.pC = zeros(n,Config.Etapas);
	Var.ClNI.pC(ClNI,:) = pCClNI;
	Var.ClNI.qC = Var.ClNI.pC * 0;
	Var.ClNI.qC(ClNI,:) = qCClNI;
	Var.ClNI.on = Var.ClNI.pC * 0;
	Var.ClNI.on(ClNI,:) = onClNI;
	Var.ClNI.start = Var.ClNI.pC * 0;
	Var.ClNI.start(ClNI,:) = stClNI;

	Var.Red.Bus.pC = Var.ClNI.pC;
	Var.Red.Bus.qC = Var.ClNI.qC;
end

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_mu);
opt(1,3) = sum(tfopt_lambda);
opt(1,4) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end