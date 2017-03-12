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

tnnLow = (1 + Data.Red.Bus.TapLow.*Data.Red.Bus.Ntr);
tnnTop = (1 + Data.Red.Bus.TapTop.*Data.Red.Bus.Ntr);

NcpCapL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow;
NcpCapT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop;

NcpvL = Data.Red.Bus.Ncp.*(Data.Red.Bus.uLow).^2;
NcpvT = Data.Red.Bus.Ncp.*(Data.Red.Bus.uTop).^2;

NcpCapLvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uLow).^2;
NcpCapTvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uTop).^2;
NcpCapLvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uTop).^2;
NcpCapTvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uLow).^2;


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
	expression tfopt_virt(Config.Etapas,1); 
	expression tfopt_conv(Config.Etapas,1); 
	expression fopt_expr; 


	%% Restricciones de generadores Solares
	% Variables de generadores solares
	variable sPv(nPv, Config.Etapas);
	variable xiPv(nPv, Config.Etapas); % module of square complex current in i
	variable cqPv(n, Config.Etapas);
	expression SPvNorm(nPv, Config.Etapas,2);

	tfopt_expr = tfopt_expr ...
		+ sum(Data.Cost.rhopPv .* pPv) ...
		+ sum(cqPv) ...
	;
	tfopt_conv = 1/(2*DistrInfo.Gama) * ...
			(norms(pPv -  DistrInfo.PV.pPv(Pv,:),2,2) ...
			+ norms(qPv - DistrInfo.PV.qPv(Pv,:),2,2));

	tfopt_virt = DistrInfo.muT(Pv,:) .* pPv + DistrInfo.lambdaT(Pv,:) .* qPv;

	cqPv >= - Data.Cost.rhomqPv .* qPv;
	cqPv >= Data.Cost.rhoMqPv .* qPv;

	pPv == (Data.Gen.Pv.pPvg(Pv,:) - (Data.Gen.Pv.cv(Pv,:).*sPv + Data.Gen.Pv.cr(Pv,:).*xiPv));

	SPvNorm(:,:,1) = pPv;
	SPvNorm(:,:,2) = qPv;
	sPv >= norms(SPvNorm,2,3);
	xiPv >= pPv.^2 + qPv.^2;

	sPv <= Data.Gen.Pv.sTop(Pv,:).*abs(sign(Data.Gen.Pv.pPvg(Pv,:)));
	xiPv <= Data.Gen.Pv.xiTop(Pv,:).*abs(sign(Data.Gen.Pv.pPvg(Pv,:)));


	fopt_expr = sum(tfopt_expr + tfopt_virt + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
Var.Gen.Pv.pPv = pPv;
Var.Gen.Pv.qPv = qPv;
Var.Gen.Pv.s = pPv*0;
Var.Gen.Pv.s(Pv,:) = sPv;
Var.Gen.Pv.xi = pPv*0;
Var.Gen.Pv.xi(Pv,:) = xiPv;

Var.Red.Bus.pG = Var.Gen.Pv.pPv;
Var.Red.Bus.qG = Var.Gen.Pv.qPv;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_virt);
opt(1,3) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end