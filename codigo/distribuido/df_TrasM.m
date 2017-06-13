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

	for i = 1: size(Config.Centr,1)
		cvx_solver_settings(Config.Centr{i,1}, Config.Centr{i,2});
	end

	cvx_precision high
	%% Declaracion de variables y expresiones
	%Variables generales de la red
	variable pGTras(nG, Config.Etapas);
	variable qGTras(nG, Config.Etapas);
	variable cqGTras(nG, Config.Etapas);

	expression tfopt_expr(Config.Etapas,1); 
	expression tfopt_virt(Config.Etapas,1); 
	expression tfopt_conv(Config.Etapas,1); 
	expression fopt_expr; 


	%% Funcion objetivo
    tfopt_conv = 1/(2*DistrInfo.Gama) * ...
            (norms(pGTras - DistrInfo.Tras.P(G,:),2,2) ...
            + norms(qGTras - DistrInfo.Tras.Q(G,:),2,2));
	tfopt_expr = ...
		sum(Data.Cost.piPTras(G,:).*pGTras,1) ...
		+ sum(cqGTras,1) ...
		;

	tfopt_virt = - DistrInfo.muT(G,:) .* pGTras ...
		- DistrInfo.lambdaT(G,:) .* qGTras;

	cqGTras >= - Data.Cost.piQmtras(G,:) .* qGTras;
	cqGTras >= Data.Cost.piQMtras(G,:) .* qGTras;

	%% Restricciones de generacion
	% Restricciones de Trasmision
	pGTras >= Data.Gen.Tras.pgLow(G,:);
	pGTras <= Data.Gen.Tras.pgTop(G,:);
	qGTras >= Data.Gen.Tras.qgLow(G,:);
	qGTras <= Data.Gen.Tras.qgTop(G,:);


	fopt_expr = sum(tfopt_expr + tfopt_virt + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT
Var.Red.Bus.PTras = zeros(n, Config.Etapas);
Var.Red.Bus.QTras = zeros(n, Config.Etapas);
Var.Red.Bus.PTras(G,:)	 = 	pGTras	;
Var.Red.Bus.QTras(G,:)	 = 	qGTras	;

Var.Red.Bus.pG = Var.Red.Bus.PTras;
Var.Red.Bus.qG = Var.Red.Bus.QTras;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_virt);
opt(1,3) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end