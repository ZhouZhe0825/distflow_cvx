function [Var, opt, status] = df_Pvm(Data, Config, DistrInfo)

%% Inicializacion
VertI = VertIMat(Data.Red.Branch.T);
VertJ = VertJMat(Data.Red.Branch.T);
OutBr = VertI';
InBr = VertJ';

n = size(VertI,2);
m = size(VertI,1);

G = find(Data.Gen.Tras.I == 1);
NnoG = setdiff((1:n),G);

Data.Red.Branch.Tup = triu(Data.Red.Branch.T);
[rowTup, colTup, ~] = find(Data.Red.Branch.Tup == 1);
Tind = find(Data.Red.Branch.T == 1);
Tupmind = sub2ind(size(Data.Red.Branch.T), rowTup, colTup);
Tdownmind = sub2ind(size(Data.Red.Branch.T), colTup, rowTup);
Tupm = zeros(size(Tupmind));
Tdownm = zeros(size(Tdownmind));
for i=1:length(Tupmind)
	Tupm(i) = find(Tind == Tupmind(i));
	Tdownm(i) = find(Tind == Tdownmind(i));
end

% Fotovoltaico
Pv = DistrInfo.Bus;
nPv = length(Pv);
NotPv = find(sign(sum(abs(Data.Gen.Pv.I),2)) == 0);

%% Modelo programacion matematica

tic
cvx_begin quiet

	for i = 1: size(Config.Centr,1)
		cvx_solver_settings(Config.Centr{i,1}, Config.Centr{i,2});
	end

	cvx_precision high
	%% Declaracion de variables y expresiones
	%Variables generales de la red
	variable pPv(nPv, Config.Etapas);
	variable qPv(nPv, Config.Etapas);

	expression tfopt_expr(Config.Etapas,1); 
	expression tfopt_mu(Config.Etapas,1); 
	expression tfopt_lambda(Config.Etapas,1); 
	expression tfopt_conv(Config.Etapas,1); 
	expression fopt_expr; 


	%% Restricciones de generadores Solares
	% Variables de generadores solares
	variable sPv(nPv, Config.Etapas);
	variable xiPv(nPv, Config.Etapas); % module of square complex current in i
	variable cqPv(nPv, Config.Etapas);
	expression SPvNorm(nPv, Config.Etapas,2);

	tfopt_expr = sum(Data.Cost.rhopPv(Pv,:) .* pPv) ...
		+ sum(cqPv) ...
	;
	tfopt_conv = 1/(2*DistrInfo.Gama) * ...
			(norms(pPv -  DistrInfo.PV.pPv(Pv,:),2,2) ...
			+ norms(qPv - DistrInfo.PV.qPv(Pv,:),2,2));

	tfopt_mu = DistrInfo.muT(Pv,:) .* pPv;
	tfopt_lambda = DistrInfo.lambdaT(Pv,:) .* qPv;

	cqPv >= - Data.Cost.rhomqPv(Pv,:) .* qPv;
	cqPv >= Data.Cost.rhoMqPv(Pv,:) .* qPv;

	pPv == (Data.Gen.Pv.pPvg(Pv,:) - (Data.Gen.Pv.cv(Pv,:).*sPv + Data.Gen.Pv.cr(Pv,:).*xiPv));

	SPvNorm(:,:,1) = pPv;
	SPvNorm(:,:,2) = qPv;
	sPv >= norms(SPvNorm,2,3);
	xiPv >= pPv.^2 + qPv.^2;

	sPv <= Data.Gen.Pv.sTop(Pv,:).*abs(sign(Data.Gen.Pv.pPvg(Pv,:)));
	xiPv <= Data.Gen.Pv.xiTop(Pv,:).*abs(sign(Data.Gen.Pv.pPvg(Pv,:)));


	fopt_expr = sum(tfopt_expr + tfopt_mu + tfopt_lambda + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
Var.Gen.Pv.pPv = zeros(n, Config.Etapas);
Var.Gen.Pv.pPv(Pv,:) = pPv;
Var.Gen.Pv.qPv = zeros(n, Config.Etapas);
Var.Gen.Pv.qPv(Pv,:) = qPv;
Var.Gen.Pv.cqPv = zeros(n, Config.Etapas);
Var.Gen.Pv.cqPv(Pv,:) = cqPv;

Var.Gen.Pv.s = zeros(n, Config.Etapas);
Var.Gen.Pv.s(Pv,:) = sPv;
Var.Gen.Pv.xi = zeros(n, Config.Etapas);
Var.Gen.Pv.xi(Pv,:) = xiPv;

Var.Red.Bus.pG = Var.Gen.Pv.pPv;
Var.Red.Bus.qG = Var.Gen.Pv.qPv;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_mu);
opt(1,3) = sum(tfopt_lambda);
opt(1,4) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end