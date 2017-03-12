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
	M = zeros(2,2,nSt,Config.Etapas);
	for i = 1:nSt
		for et = 1: Config.Etapas
			M(:,:,i,et) = [Data.St.Bat.m1(St(i),et) -Data.St.Bat.m2(St(i),et)/2; ...
				-Data.St.Bat.m2(St(i),et)/2 Data.St.Bat.m1(St(i),et)];
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
	expression tfopt_virt(Config.Etapas,1); 
	expression tfopt_conv(Config.Etapas,1); 
	expression fopt_expr; 



	%% Restricciones de almacenamiento de bateria
	variable sStb(nSt, Config.Etapas);
	variable pStgb(nSt, Config.Etapas);
	variable xiStb(nSt, Config.Etapas);
	variable EStb(nSt, Config.Etapas);
	variable DlEStb(nSt, Config.Etapas);

	expression StbNorm(nSt, Config.Etapas,2);
	expression EStbAnt(nSt, Config.Etapas);
	expression cStb(n, Config.Etapas);

	DlEStb <= 0;
	DlEStb <= EStb - Data.St.Bat.ETop*Data.St.Bat.kapa;


	cStb = (Data.St.Bat.wOm(St,:) + Data.St.Bat.m3(St,:).*(DlEStb.^2)); % falta termino de m2

	% Modelado de Config.Etapas
	for et = 1: Config.Etapas
		for i = 1:nSt
			j = St(i);
			if et == 1
				cStb(j,et) = cStb(j,et) + [pStgb(i,et) 0] * M(:,:,i,et) * [pStgb(i,et); 0];
			else
				cStb(j,et) = cStb(j,et) + [pStgb(i,et) pStgb(i,et-1)] * M(:,:,i,et) * [pStgb(i,et); pStgb(i,et-1)];
			end
		end
	end

	tfopt_expr = tfopt_expr + sum(cStb,1);

	tfopt_expr = ...
			sum(Data.St.Bat.I(:,Config.Etapas).*Data.St.Bat.beta(:,Config.Etapas).* ...
			((Data.St.Bat.ETop(:,Config.Etapas) - EStb(:,Config.Etapas)*Data.St.Bat.gama).^2) ...
			+ Data.St.Bat.wU(:,Config.Etapas),1)./Config.Etapas;

	tfopt_virt = DistrInfo.muT(St,:) .* pStb + DistrInfo.lambdaT(St,:) .* qStb;

	tfopt_conv = 1/(2*DistrInfo.Gama) * ...
			(norms(pStb - DistrInfo.Alm.pStb(St,:),2,2) ...
			+ norms(qStb - DistrInfo.Alm.qStb(St,:),2,2));

	fopt_expr = sum(tfopt_expr) + sum(tfopt_virt) + sum(tfopt_conv);

	EStbAnt(:,1) = Data.St.Bat.EIni(St,1);
	EStbAnt(:,(2:Config.Etapas)) = EStb(:,(1:Config.Etapas-1));

	pStb == pStgb - (Data.St.Bat.cv(St,:).*sStb + Data.St.Bat.cr(St,:).*xiStb);
	EStb == (1-Data.St.Bat.epsilon(St,:)).*EStbAnt - Data.St.Bat.eta(St,:).*pStgb*Data.dt;

	StbNorm(:,:,1) = pStgb;
	StbNorm(:,:,2) = qStb;
	sStb >= norms(StbNorm,2,3);

	xiStb >= pStgb.^2 + qStb.^2;

	pStb <= Data.St.Bat.pgTop(St,:);
	sStb <= Data.St.Bat.sTop(St,:);
	xiStb <= Data.St.Bat.xiTop(St,:);
	EStb <= Data.St.Bat.ETop(St,:);

	pStb >= Data.St.Bat.pgLow(St,:);
	EStb >= Data.St.Bat.ELow(St,:);

	fopt_expr = sum(tfopt_expr + tfopt_virt + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT

% Baterias
Var.St.Bat.pStb = zeros(n, Config.ETapas);
Var.St.Bat.pStb(St,:) = pStb;
Var.St.Bat.pStgb = pStb*0;
Var.St.Bat.pStgb(St,:) = pStgb;
Var.St.Bat.qStb = pStb*0;
Var.St.Bat.qStb(St,:) = qStb;
Var.St.Bat.sStb = pStb*0;
Var.St.Bat.sStb(St,:) = sStb;
Var.St.Bat.xiStb = pStb*0;
Var.St.Bat.xiStb(St,:) = xiStb;
Var.St.Bat.EStb = pStb*0;
Var.St.Bat.EStb(St,:) = EStb;


Var.Red.Bus.pG = Var.St.Bat.pStb;
Var.Red.Bus.qG = Var.St.Bat.qStb;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_virt);
opt(1,3) = sum(tfopt_conv);


status = cvx_status;
cvx_clear
end