function [Var, opt, status] = df_TrasM(Data, Config, DistrInfo)

%% Inicializacion
VertI = VertIMat(Data.Red.Branch.T);
VertJ = VertJMat(Data.Red.Branch.T);
OutBr = VertI';
InBr = VertJ';

n = size(VertI,2);
m = size(VertI,1);

G = DistrInfo.Bus;
nG = length(G);


%% Modelo programacion matematica
tic
cvx_begin quiet

	for i = 1: size(Config.SubP,1)
		cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
	end

	cvx_precision high
	%% Declaracion de variables y expresiones
	%Variables generales de la red
	variable pGTras(nG, Config.Etapas);
	variable qGTras(nG, Config.Etapas);
	variable cqGTras(nG, Config.Etapas);

	expression tfopt_expr(Config.Etapas,1); 
	expression tfopt_mu(Config.Etapas,1); 
	expression tfopt_lambda(Config.Etapas,1); 
	expression tfopt_conv(Config.Etapas,1); 
	expression fopt_expr; 


	%% Funcion objetivo
	tfopt_conv = 1/(2*DistrInfo.Gama) * ...
			(norms(pGTras - DistrInfo.Tras.P(G,:),2,1) ...
			+ norms(qGTras - DistrInfo.Tras.Q(G,:),2,1));
	tfopt_expr = ...
		sum(Data.Cost.piPTras(G,:).*pGTras,1) ...
		+ sum(cqGTras,1) ...
		;

	tfopt_mu = sum(- DistrInfo.muT(G,:) .* pGTras,1);
	tfopt_lambda = sum(- DistrInfo.lambdaT(G,:) .* qGTras,1);

	cqGTras >= - Data.Cost.piQmtras(G,:) .* qGTras;
	cqGTras >= Data.Cost.piQMtras(G,:) .* qGTras;

	%% Restricciones de generacion
	% Restricciones de Trasmision
	pGTras >= Data.Gen.Tras.pgLow(G,:);
	pGTras <= Data.Gen.Tras.pgTop(G,:);
	qGTras >= Data.Gen.Tras.qgLow(G,:);
	qGTras <= Data.Gen.Tras.qgTop(G,:);


	fopt_expr = sum(tfopt_expr + tfopt_mu + tfopt_lambda + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT
Var.Red.Bus.PTras = zeros(n, Config.Etapas);
Var.Red.Bus.QTras = zeros(n, Config.Etapas);
Var.Red.Bus.cQTras = zeros(n, Config.Etapas);
Var.Red.Bus.PTras(G,:)	 = 	pGTras	;
Var.Red.Bus.QTras(G,:)	 = 	qGTras	;
Var.Red.Bus.cQTras(G,:)  =  cqGTras ;


Var.Red.Bus.pG = Var.Red.Bus.PTras;
Var.Red.Bus.qG = Var.Red.Bus.QTras;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_mu);
opt(1,3) = sum(tfopt_lambda);
opt(1,4) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end