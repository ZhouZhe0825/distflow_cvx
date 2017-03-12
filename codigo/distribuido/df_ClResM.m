function [Var, opt, status] = df_ClResM(Data, Config, DistrInfo)

%% Inicializacion
VertI = VertIMat(Data.Red.Branch.T);
VertJ = VertJMat(Data.Red.Branch.T);
OutBr = VertI';
InBr = VertJ';

n = size(VertI,2);
m = size(VertI,1);

DistrInfoM.ClRes.v	 = 	replicateMat3DApp(	DistrInfo.ClRes.v	,2);
DistrInfoM.lambdaT	 = 	replicateMat3DApp(	DistrInfo.lambdaT	,2);
DistrInfoM.muT	 = 	replicateMat3DApp(	DistrInfo.muT	,2);


% Aire Acondicionado
AC = find(Data.St.AC.I == 1);
nAC = length(AC);

%% Modelo programacion matematica

tic
cvx_begin quiet

	for i = 1: size(Config.SubP,1)
		cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
	end

	cvx_precision high
	%% Declaracion de variables y expresiones
	%Variables generales de la red
	variable pCApp(n, Config.Etapas, 2); % real power demand in i
	variable qCApp(n, Config.Etapas, 2); % real power demand in i
	variable pCn(n, Config.Etapas, 2); % real power demand in i
	variable pCClRes(n, Config.Etapas); % real power demand in i
	variable qCClRes(n, Config.Etapas); % real power demand in i

	expression vApp(n, Config.Etapas, 2);

	expression tfopt_expr(Config.Etapas,1); 
	expression tfopt_virt(Config.Etapas,1); 
	expression tfopt_conv(Config.Etapas,1); 
	expression fopt_expr; 


	%% Funcion objetivo
	tfopt_expr = ...
		sum(Data.Util.betaT(:,:,1).*((pCn(:,:,1) - Data.Util.pzCnPref(:,:,1)).^2),1) ...
		;

	tfopt_virt = sum(sum(DistrInfoM.muT(DistrInfo.Bus,:,:) .* pCApp(DistrInfo.Bus,:,:) + DistrInfoM.lambdaT(DistrInfo.Bus,:,:) .* qCApp(DistrInfo.Bus,:,:),3),1);

	tfopt_conv = 1/(2*DistrInfo.Gama) * ...
		sum(sum(norms(pCApp(DistrInfo.Bus,:,:) - DistrInfo.ClRes.pCApp(DistrInfo.Bus,:,:),2,2) ...
		+ norms(qCApp(DistrInfo.Bus,:,:) - DistrInfo.ClRes.qCApp(DistrInfo.Bus,:,:),2,2),3),1);

	%% Restricciones de clientes residenciales
	uLowApp = zeros(size(Data.Red.Bus.uLow,1),size(Data.Red.Bus.uLow,2), 2);
	uTopApp = uLowApp;

	for app = 1:2
		vApp(:,:,app) = DistrInfo.ClRes.v;
		uLowApp(:,:,app) = Data.Red.Bus.uLow(:,:);
		uTopApp(:,:,app) = Data.Red.Bus.uTop(:,:);
	end

	pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnLow.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn; 
	pCApp >= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnTop.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

	pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnTop.*vApp + pCn.*uLowApp.^2 - Data.Util.pzCnTop.*uLowApp.^2) + (1-Data.Red.Bus.alpha).* pCn;
	pCApp <= Data.Red.Bus.alpha.*(Data.Util.pzCnLow.*vApp + pCn.*uTopApp.^2 - Data.Util.pzCnLow.*uTopApp.^2) + (1-Data.Red.Bus.alpha).* pCn;

	qCApp == pCApp.*Data.Util.tgPhi;

	pCn >= Data.Util.pzCnLow;
	pCn <= Data.Util.pzCnTop;

	pCClRes >= sum(pCApp,3);

	qCClRes >= sum(qCApp,3);

	%% Restricciones de Aire Acondicionado
	if nAC > 0 
		% Temperatura Aire Acondicionado
		variable Tvar(n,Config.Etapas);
		expression TvarAnt(n,Config.Etapas);
		expression pCAC(n,Config.Etapas);

		TvarAnt(:,1) = Data.St.AC.tempIni;
		TvarAnt(:,(2:Config.Etapas)) = Tvar(:,(1:Config.Etapas-1));
		pCAC(:,:) = 0;
		pCAC(AC,:) = pCApp(AC,:,2);

		Tvar(AC,:) - TvarAnt(AC,:) + Data.St.AC.epsilon(AC,:).*(TvarAnt(AC,:) - Data.temp(AC,:))*Data.dt - Data.St.AC.eta(AC,:).*pCAC(AC,:)*Data.dt == 0;

		Tvar(AC,:) >= Data.St.AC.tempLow(AC,:);
		Tvar(AC,:) <= Data.St.AC.tempTop(AC,:);
		tfopt_expr = tfopt_expr + sum(Data.St.AC.beta.*(Tvar - Data.St.AC.tempPref).^2,1);

	end

	fopt_expr = sum(tfopt_expr + tfopt_virt + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT

Var.ClRes.pCApp = full(pCApp);
Var.ClRes.qCApp = full(qCApp);
Var.ClRes.pC = full(pCClRes);
Var.ClRes.qC = full(qCClRes);

% Aire Acondicionado
if nAC > 0 
	Var.ClRes.Tvar = Tvar;
end

Var.Red.Bus.pC = Var.ClRes.pC;
Var.Red.Bus.qC = Var.ClRes.qC;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_virt);
opt(1,3) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end

function [A] = replicateMat3DApp(V, App)
	A = zeros(size(V,1), size(V,2), App);
	for app = 1:App
		A(:,:,app) = V;
	end
end
