function [Var_dist_conE, Var_centr, Var_F, opt_dist_conE, opt_centr, opt_F, status, Data_d, Ev] = llamarDistribuidoM(Data, Config, Var_centr, opt_centr, Var_ini)

	DataM = Data;
	[DataM] = reshapeDataM(DataM, Config);

	if isempty(Var_centr) && isempty(Var_ini)

		ConfigIni = Config;
		ConfigIni.Centr = Config.Ini;

		leyenda = 'Primera vuelta centralizado'
		[Var_centr, opt_centr, status] = distflowCentralizadoM(DataM, Config);
		leyenda = ['Centralizado - ' status]

		if isNotSolved(status)
			save(Config.workspace_var_file);
			return;
		end

		[Var_ini, opt_ini, status] = distflowCentralizadoM(DataM, ConfigIni);
		if isNotSolved(status)
			save(Config.workspace_var_file);
			return;
		end
		leyenda = ['Inicial - ' status]

		DataM_f = DataM;

		leyenda = 'Segunda vuelta centralizado'


%		 DataM_f.Fixed.Cap = Var_ini.Red.Bus.Cap;
%		 DataM_f.Fixed.Tap = Var_ini.Red.Bus.Tap; 
%		 DataM_f.Fixed.stCh = Var_ini.ClNI.start;
%		 DataM_f.Fixed.onCh = Var_ini.ClNI.on;
%		 DataM_f.Fixed.y = Var_ini.Red.Branch.y;
%		 DataM_f.Fixed.z = Var_ini.Red.Branch.z;

		DataM_f.Fixed.Cap = round(Var_ini.Red.Bus.Cap);
		DataM_f.Fixed.Tap = round(Var_ini.Red.Bus.Tap);
        if isfield(Var_ini, 'ClNI')
            DataM_f.Fixed.stCh = round(Var_ini.ClNI.start);
            DataM_f.Fixed.onCh = round(Var_ini.ClNI.on);
        end
		DataM_f.Fixed.y = round(Var_ini.Red.Branch.y);
		DataM_f.Fixed.z = round(Var_ini.Red.Branch.z);
        
		[Var_ini_, opt_ini, status] = distflowCentralizadoM(DataM_f, ConfigIni);

% 		difCap_ini = max(max(abs(DataM_f.Fixed.Cap - Var_ini.Red.Bus.Cap)));
% 		difTap_ini = max(max(abs(DataM_f.Fixed.Tap - Var_ini.Red.Bus.Tap)));
% 		difstCh_ini = max(max(abs(DataM_f.Fixed.stCh - Var_ini.ClNI.start)));
% 		difonCh_ini = max(max(abs(DataM_f.Fixed.onCh - Var_ini.ClNI.on)));
% 		dify_ini = max(max(abs(DataM_f.Fixed.y - Var_ini.Red.Branch.y)));
% 		difz_ini = max(max(abs(DataM_f.Fixed.z - Var_ini.Red.Branch.z)));

		leyenda = ['Inicial Duales - ' status]

		leyenda = 'Fin centralizado'

		Var_ini.Dual.dPn = Var_ini_.Dual.dPn;
		Var_ini.Dual.dQn = Var_ini_.Dual.dQn;

		[Var_ini] = randStart(Var_ini);
		save(Config.workspace_var_file);
	end
	tdistr = tic;
	nodos = size(Data.Red.Branch.T,1);
	DistrInfo.Bus = [1:nodos];
	DistrInfo.Gama = Config.Distr.gama_ini;

    TrasNodes = find(Data.Gen.Tras.I == 1);
	DfigNodes = find(Data.Gen.DFIG.I == 1);
	PvNodes = find(Data.Gen.Pv.I == 1);
	AlmNodes = find(Data.St.Bat.I == 1);
	ClResNodes = Data.Red.Bus.indCons;
	ClNINodes = find(Data.ClNI.I == 1);

	ConsGenNodes = [TrasNodes];
	ConsGenNodes = union(ConsGenNodes, DfigNodes);
	ConsGenNodes = union(ConsGenNodes, PvNodes);
	ConsGenNodes = union(ConsGenNodes, AlmNodes);
	ConsGenNodes = union(ConsGenNodes, ClResNodes);
	ConsGenNodes = union(ConsGenNodes, ClNINodes);

	OpSisNodes = setdiff((1:size(Data.Red.Branch.T,1)), ConsGenNodes);
	DataM.Red.Bus.indCG0 = OpSisNodes;



	[DistrInfo] = actualizarDistrInfo(Var_ini, DataM, TrasNodes, DfigNodes, ClNINodes, PvNodes, AlmNodes, DistrInfo);

	% DistrInfo general

	DistrInfo.muT = zeros(size(Var_ini.Dual.dPn,1),1,size(Var_ini.Dual.dPn,2));
	DistrInfo.lambdaT = DistrInfo.muT;
	DistrInfo.muT = Var_ini.Dual.dPn;
	DistrInfo.lambdaT = Var_ini.Dual.dQn;

	DistrInfo.mu = DistrInfo.muT;
	DistrInfo.lambda = DistrInfo.lambdaT;

	Var_dist_conE = Var_ini;

	Ev.opt_Tras	 = zeros(Config.Distr.TopIT*2+1+1,3);
	Ev.opt_Dfig	 = zeros(Config.Distr.TopIT*2+1+1,3);
	Ev.opt_Pv	 = zeros(Config.Distr.TopIT*2+1+1,3);
	Ev.opt_ClNI	 = zeros(Config.Distr.TopIT*2+1+1,3);
	Ev.opt_ClRes	 = zeros(Config.Distr.TopIT*2+1+1,3);
	Ev.opt_Alm	 = zeros(Config.Distr.TopIT*2+1+1,3);
	Ev.opt_OpSis	 = zeros(Config.Distr.TopIT*2+1+1,3);
	Ev.opt	 = zeros(Config.Distr.TopIT*2+1,3);
	Ev.DifP	 = zeros(nodos,Config.Etapas,Config.Distr.TopIT*2+1);
	Ev.DifQ	 = zeros(nodos,Config.Etapas,Config.Distr.TopIT*2+1);
	Ev.errP	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;
	Ev.errQ	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;
	Ev.mu	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * 0;
	Ev.lambda	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * 0;
	Ev.pN	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;
	Ev.qN	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;
	Ev.pG	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;
	Ev.qG	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;
	Ev.pC	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;
	Ev.qC	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;
	Ev.errMu	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;
	Ev.errLambda	 = ones(nodos, Config.Etapas, Config.Distr.TopIT*2+1) * Inf;


	Data_c = DataM;
	Data_d = DataM;

	wbSol = waitbar(0,'1','Name',['Calculo de solucion distr']);
	[Var_dist_conE, opt_dist_conE, status, it, DistrInfo, Ev, wbSol] = DistrLoop(0, false, DistrInfo, Config, DataM, Ev, wbSol);
	if isNotSolved(status)
		save(Config.workspace_var_file);
		return;
    end
    
    it2 = 0;    
% % 	save(Config.workspace_var_file);
% % 
% % 	DataM_d = DataM;
% % 
% % 	DataM_d.Fixed.Cap = round(Var_dist_conE.Red.Bus.Cap);
% % 	DataM_d.Fixed.Tap = round(Var_dist_conE.Red.Bus.Tap);
% % 	if isfield(Var_dist_conE, 'ClNI')
% % 		DataM_d.Fixed.stCh = round(Var_dist_conE.ClNI.start);
% % 		DataM_d.Fixed.onCh = round(Var_dist_conE.ClNI.on);
% % 	end
% % 	DataM_d.Fixed.y = round(Var_dist_conE.Red.Branch.y);
% % 	DataM_d.Fixed.z = round(Var_dist_conE.Red.Branch.z);
% % 
% % % 	Ev.difCap = max(max(abs(DataM_d.Fixed.Cap - Var_dist_conE.Red.Bus.Cap)));
% % % 	Ev.difTap = max(max(abs(DataM_d.Fixed.Tap - Var_dist_conE.Red.Bus.Tap)));
% % % 	if isfield(Var_dist_conE, 'ClNI')
% % % 		Ev.difstCh = max(max(abs(DataM_d.Fixed.stCh - Var_dist_conE.ClNI.start)));
% % % 		Ev.difonCh = max(max(abs(DataM_d.Fixed.onCh - Var_dist_conE.ClNI.on)));
% % % 	end
% % % 	Ev.dify = max(max(abs(DataM_d.Fixed.y - Var_dist_conE.Red.Branch.y)));
% % % 	Ev.difz = max(max(abs(DataM_d.Fixed.z - Var_dist_conE.Red.Branch.z)));
% % 
% % 	
% % 	[Var_F, opt_F, status, it2, DistrInfo, Ev, wbSol] = DistrLoop(it, true, DistrInfo, Config, DataM_d, Ev, wbSol);
% % 	if isNotSolved(status)
% % 		save(Config.workspace_var_file);
% % 		return;
% % 	end
% % 

    Var_F = Var_dist_conE;
    opt_F = opt_dist_conE;
    
    delete(wbSol);
	Ev.opt = Ev.opt((1:it2+it-1),:);
	Ev.DifP = Ev.DifP(:,:,(1:it2+it-1));
	Ev.DifQ = Ev.DifQ(:,:,(1:it2+it-1));

	Ev.mu = Ev.mu(:,:,(1:it2+it-1));
	Ev.lambda = Ev.lambda(:,:,(1:it2+it-1));

	save(Config.workspace_var_file);
	toc(tdistr)
end

function [Var_mod] = randStart(Var)
	Var_mod = Var;

	amplitudD = .15;
	amplitud = .1;
%	 
%	 rnd = 1 + amplitud * (rand - .5);
	rnd = 1;


	if isfield(Var_mod,'Gen')
		if isfield(Var_mod.Gen, 'Dfig')
			Var_mod.Gen.Dfig.pWi = ones(size(Var.Gen.Dfig.pWi)) * rnd;
			Var_mod.Gen.Dfig.qWi = ones(size(Var.Gen.Dfig.qWi)) * rnd;

			Var_mod.Gen.Dfig.pWi = Var_mod.Gen.Dfig.pWi .* Var.Gen.Dfig.pWi;
			Var_mod.Gen.Dfig.qWi = Var_mod.Gen.Dfig.qWi .* Var.Gen.Dfig.qWi;
		end
		if isfield(Var_mod.Gen, 'Pv')
			Var_mod.Gen.Pv.pPv = ones(size(Var.Gen.Pv.pPv)) * rnd;
			Var_mod.Gen.Pv.qPv = ones(size(Var.Gen.Pv.qPv)) * rnd;

			Var_mod.Gen.Pv.pPv = Var_mod.Gen.Pv.pPv .* Var.Gen.Pv.pPv;
			Var_mod.Gen.Pv.qPv = Var_mod.Gen.Pv.qPv .* Var.Gen.Pv.qPv;
		end
	end


	if isfield(Var_mod,'St')
		if isfield(Var_mod.Gen, 'Bat')
			Var_mod.St.Bat.pStb = ones(size(Var.St.Bat.pStb)) * rnd;
			Var_mod.St.Bat.qStb = ones(size(Var.St.Bat.qStb)) * rnd;

			Var_mod.St.Bat.pStb = Var_mod.St.Bat.pStb .* Var.St.Bat.pStb;
			Var_mod.St.Bat.qStb = Var_mod.St.Bat.qStb .* Var.St.Bat.qStb;
		end
	end

	Var_mod.Red.Bus.PTras = ones(size(Var.Red.Bus.PTras)) * rnd;
	Var_mod.Red.Bus.QTras = ones(size(Var.Red.Bus.QTras)) * rnd;
	Var_mod.Red.Bus.qCp = ones(size(Var.Red.Bus.qCp)) * rnd;
	Var_mod.Red.Bus.v = ones(size(Var.Red.Bus.v)) * rnd;
	Var_mod.Red.Bus.pN = ones(size(Var.Red.Bus.pN)) * rnd;
	Var_mod.Red.Bus.qN = ones(size(Var.Red.Bus.qN)) * rnd;

	Var_mod.ClRes.pC = ones(size(Var.ClRes.pC)) * rnd;
	Var_mod.ClRes.qC = ones(size(Var.ClRes.qC)) * rnd;
	Var_mod.ClRes.pCApp = ones(size(Var.ClRes.pCApp)) * rnd;
	Var_mod.ClRes.qCApp = ones(size(Var.ClRes.qCApp)) * rnd;

	Var_mod.Red.Branch.P = ones(size(Var.Red.Branch.P)) * rnd;
	Var_mod.Red.Branch.Q = ones(size(Var.Red.Branch.Q)) * rnd;


	Var_mod.Red.Bus.PTras = Var_mod.Red.Bus.PTras .* Var.Red.Bus.PTras;
	Var_mod.Red.Bus.QTras = Var_mod.Red.Bus.QTras .* Var.Red.Bus.QTras;
	Var_mod.Red.Bus.qCp = Var_mod.Red.Bus.qCp .* Var.Red.Bus.qCp;
	Var_mod.Red.Bus.v = Var_mod.Red.Bus.v .* Var.Red.Bus.v;
	Var_mod.Red.Bus.pN = Var_mod.Red.Bus.pN .* Var.Red.Bus.pN;
	Var_mod.Red.Bus.qN = Var_mod.Red.Bus.qN .* Var.Red.Bus.qN;

	Var_mod.ClRes.pC = Var_mod.ClRes.pC .* Var.ClRes.pC;
	Var_mod.ClRes.qC = Var_mod.ClRes.qC .* Var.ClRes.qC;
	Var_mod.ClRes.pCApp = Var_mod.ClRes.pCApp .* Var.ClRes.pCApp;
	Var_mod.ClRes.qCApp = Var_mod.ClRes.qCApp .* Var.ClRes.qCApp;

	Var_mod.Red.Branch.P = Var_mod.Red.Branch.P .* Var.Red.Branch.P;
	Var_mod.Red.Branch.Q = Var_mod.Red.Branch.Q .* Var.Red.Branch.Q;

	Var_mod.Dual.dPn = (1 + amplitudD * (rand(size(Var.Dual.dPn)) - .5)).*Var.Dual.dPn;
	Var_mod.Dual.dQn = (1 + amplitudD * (rand(size(Var.Dual.dPn)) - .5)).*Var.Dual.dQn;

end

function [DistrInfo] = actualizarDistrInfo(Var, Data, TrasNodes, DfigNodes, ClNINodes, PvNodes, AlmNodes, DistrInfo)

	% DistrInfo Administrador
	DistrInfo.Adm.pG = Var.Red.Bus.pG;
	DistrInfo.Adm.qG = Var.Red.Bus.qG;
	DistrInfo.Adm.pN = Var.Red.Bus.pN;
	DistrInfo.Adm.qN = Var.Red.Bus.qN;
	DistrInfo.Adm.pC = Var.Red.Bus.pC;
	DistrInfo.Adm.qC = Var.Red.Bus.qC;

	if ~isempty(DfigNodes)
		% DistrInfo Dfig
		DistrInfo.Dfig.v = Var.Red.Bus.v;
		DistrInfo.Dfig.pWi = Var.Gen.Dfig.pWi;
		DistrInfo.Dfig.qWi = Var.Gen.Dfig.qWi;
		DistrInfo.Dfig.ind = Data.Gen.DFIG.I;
		DistrInfo.Dfig.ind(DfigNodes) = (1:length(DfigNodes))';
	end

	if ~isempty(PvNodes)
		% DistrInfo PV
		DistrInfo.PV.pPv = Var.Gen.Pv.pPv;
		DistrInfo.PV.qPv = Var.Gen.Pv.qPv;
	end

	if ~isempty(AlmNodes);
		% DistrInfo Almacenamiento
		DistrInfo.Alm.pStb = Var.St.Bat.pStb;
		DistrInfo.Alm.qStb = Var.St.Bat.qStb;
	end
	if ~isempty(ClNINodes)
		% DistrInfo Clientes No Interrumpibles
		DistrInfo.ClNI.v = Var.Red.Bus.v;
		DistrInfo.ClNI.pC = Var.ClNI.pC;
		DistrInfo.ClNI.qC = Var.ClNI.qC;
		DistrInfo.ClNI.on = Var.ClNI.on;
		DistrInfo.ClNI.start = Var.ClNI.start;
	end

	% DistrInfo Trasmision
	DistrInfo.Tras.P = Var.Red.Bus.PTras;
	DistrInfo.Tras.Q = Var.Red.Bus.QTras;

	% DistrInfo Clientes Residenciales
	DistrInfo.ClRes.v = Var.Red.Bus.v;
	DistrInfo.ClRes.pCApp = Var.ClRes.pCApp;
	DistrInfo.ClRes.qCApp = Var.ClRes.qCApp;

	% DistrInfo Operador de la red
	DistrInfo.OpSis.pN = Var.Red.Bus.pN;
	DistrInfo.OpSis.qN = Var.Red.Bus.qN;
	DistrInfo.OpSis.qCp = Var.Red.Bus.qCp;
	DistrInfo.OpSis.P = Var.Red.Branch.P;
	DistrInfo.OpSis.Q = Var.Red.Branch.Q;

end

% function [conv, stepGama, pNEv, qNEv, sNEv] = converg(pN, qN, sN, pNEv, qNEv, sNEv, optEv, DistrInfo, step, it, tolConv, tolDual)
function [conv, Ev] = converg(Ev, it, tolConv)
	if it == 1
		Ev.errMu(:,:,it) = Inf;
		Ev.errLambda(:,:,it) = Inf;
		conv = false;
	else
		ant = it - 1;
		Ev.errMu(:,:,it) = abs(Ev.mu(:,:,ant) - Ev.mu(:,:,it)) ./ abs(Ev.mu(:,:,it));
		Ev.errLambda(:,:,it) = abs(Ev.lambda(:,:,ant) - Ev.lambda(:,:,it)) ./ abs(Ev.lambda(:,:,it));
		conv = max(max(Ev.errMu(:,:,it))) <= tolConv && max(max(Ev.errLambda(:,:,it))) <= tolConv;
	end
end

function [A] = expandMat(V)
	A = zeros(size(V,1), 1, size(V,2));
	A(:,1,:) = V;
end

function [B] = contractMat(V)
	B = zeros(size(V,1), size(V,3));
	B(:,:) = V(:,1,:);
end

function [A] = expandMatApp(V)
	A = zeros(size(V,1), 1, size(V,2), size(V,3));
	A(:,1,:,:) = V;
end

function [Var_curr, opt_curr, status, it, DistrInfo, Ev, wbSol] = DistrLoop(ini_it, Fixed, DistrInfo, Config, Data, Ev, wbSol)
	conv = false;
	it = 1;
    TrasNodes = find(Data.Gen.Tras.I == 1);
	DfigNodes = find(Data.Gen.DFIG.I == 1);
	PvNodes = find(Data.Gen.Pv.I == 1);
	ClResNodes = Data.Red.Bus.indCons;
	AlmNodes = find(Data.St.Bat.I == 1);
	OpSisNodes = (1:size(Data.Red.Branch.T,1));
	ClNINodes = find(Data.ClNI.I == 1);

	while ~conv && it <= Config.Distr.TopIT

		it
		print_it = it + ini_it;

		tic
		opt_curr = zeros(1,3);
		% 1- Actualizacion del muT y lambdaT
		DistrInfo.muT = DistrInfo.mu + DistrInfo.Gama * (DistrInfo.Adm.pC - DistrInfo.Adm.pG - DistrInfo.Adm.pN);
		DistrInfo.lambdaT = DistrInfo.lambda + DistrInfo.Gama * (DistrInfo.Adm.qC - DistrInfo.Adm.qG - DistrInfo.Adm.qN);

		% 2- Invocacion a todos los participantes con lambda y mu
		[Var_curr, opt_curr, Ev.opt_Tras, Ev.opt_Dfig, Ev.opt_Pv, Ev.opt_ClNI, Ev.opt_ClRes, Ev.opt_Alm, Ev.opt_OpSis, status] = ...
			solve_all_subproblems(Data, Config, DistrInfo, print_it, wbSol, Ev.opt_Tras, Ev.opt_Dfig, Ev.opt_Pv, Ev.opt_ClNI, Ev.opt_ClRes, Ev.opt_Alm, Ev.opt_OpSis, opt_curr, Fixed);

		% 3-Segundo Construccion de lambda y mu con los consumos,
		% generadores, transmision y capacitores
		[DistrInfo] = actualizarDistrInfo(Var_curr, Data, TrasNodes, DfigNodes, ClNINodes, PvNodes, AlmNodes, DistrInfo);
		% gamaEv(print_it) = DistrInfo.Gama;

		Ev.difP(:,:,print_it) = DistrInfo.Adm.pC - DistrInfo.Adm.pG - DistrInfo.Adm.pN;
		Ev.difQ(:,:,print_it) = DistrInfo.Adm.qC - DistrInfo.Adm.qG - DistrInfo.Adm.qN;

		DistrInfo.mu = DistrInfo.mu + DistrInfo.Gama * (DistrInfo.Adm.pC - DistrInfo.Adm.pG - DistrInfo.Adm.pN);
		DistrInfo.lambda = DistrInfo.lambda + DistrInfo.Gama * (DistrInfo.Adm.qC - DistrInfo.Adm.qG - DistrInfo.Adm.qN);

		Ev.pN(:,:,print_it) = DistrInfo.Adm.pN;
		Ev.qN(:,:,print_it) = DistrInfo.Adm.qN;
		Ev.pG(:,:,print_it) = DistrInfo.Adm.pG;
		Ev.qG(:,:,print_it) = DistrInfo.Adm.qG;
		Ev.pC(:,:,print_it) = DistrInfo.Adm.pC;
		Ev.qC(:,:,print_it) = DistrInfo.Adm.qC;

		Ev.mu(:,:,print_it) = DistrInfo.mu;
		Ev.lambda(:,:,print_it) = DistrInfo.lambda;

		Ev.errP(:,:,print_it) = abs(Ev.pC(:,:,print_it) - Ev.pG(:,:,print_it)) ./ abs(Ev.pN(:,:,print_it));
		Ev.errQ(:,:,print_it) = abs(Ev.qC(:,:,print_it) - Ev.qG(:,:,print_it)) ./ abs(Ev.qN(:,:,print_it));

		Ev.opt(print_it,:) = opt_curr;

		[conv, Ev] = converg(Ev, print_it, Config.Distr.tol);
% 		save([Config.workspace_var_file '_' num2str(print_it)]);
		if isNotSolved(status)
			break;
		end

		it = it + 1;
		toc
	end
end

function [Var, opt, opt_Tras, opt_Dfig, opt_Pv, opt_ClNI, opt_ClRes, opt_Alm, opt_OpSis, status] = solve_all_subproblems(Data, Config, DistrInfo, it, wbSol, opt_Tras, opt_Dfig, opt_Pv, opt_ClNI, opt_ClRes, opt_Alm, opt_OpSis, opt, Fixed)

    TrasNodes = find(matOverTime(Data.Gen.Tras.I) == 1);
	DfigNodes = find(matOverTime(Data.Gen.DFIG.I) == 1);
	PvNodes = find(matOverTime(Data.Gen.Pv.I) == 1);
	AlmNodes = find(matOverTime(Data.St.Bat.I) == 1);
	ClResNodes = Data.Red.Bus.indCons;
	ClNINodes = find(matOverTime(Data.ClNI.I) == 1);

	OpSisNodes = (1:size(Data.Red.Branch.T,1));

	totalProb = 2 + length(DfigNodes) + length(PvNodes) + length(ClResNodes) + length(ClNINodes) + length(AlmNodes);
	% Trasmision
	previo = 0;
	[Var, opt_Tras, opt, status_Tras] = solve_problem('Tras', previo, totalProb, opt_Tras, opt, @df_TrasM, it, TrasNodes, true, wbSol, [], Data, Config, DistrInfo, Fixed);

	% Dfig
	previo = previo+1;
	[Var, opt_Dfig, opt, status_Dfig] = solve_problem('Dfig', previo, totalProb, opt_Dfig, opt, @df_DfigM, it, DfigNodes, false, wbSol, Var, Data, Config, DistrInfo, Fixed);

	% Pv
	previo = previo+length(DfigNodes);
	[Var, opt_Pv, opt, status_Pv] = solve_problem('Pv', previo, totalProb, opt_Pv, opt, @df_PvM, it, PvNodes, false, wbSol, Var, Data, Config, DistrInfo, Fixed);

	% ClNI
	previo = previo+length(PvNodes);
	[Var, opt_ClNI, opt, status_ClNI] = solve_problem('ClNI', previo, totalProb, opt_ClNI, opt, @df_ClNIM, it, ClNINodes, false, wbSol, Var, Data, Config, DistrInfo, Fixed);

	% Clientes residenciales
	previo = previo+length(ClNINodes);
	[Var, opt_ClRes, opt, status_ClRes] = solve_problem('ClRes', previo, totalProb, opt_ClRes, opt, @df_ClResM, it, ClResNodes, true, wbSol, Var, Data, Config, DistrInfo, Fixed);

	% Almacenamiento de baterias
	previo = previo+length(ClResNodes);
	[Var, opt_Alm, opt, status_Alm] = solve_problem('Alm', previo, totalProb, opt_Alm, opt, @df_AlmM, it, AlmNodes, false, wbSol, Var, Data, Config, DistrInfo, Fixed);

	% Red
	previo = previo+length(AlmNodes);
	[Var, opt_OpSis, opt, status_OpSis] = solve_problem('OpSis', previo, totalProb, opt_OpSis, opt, @df_OpSisM, it, OpSisNodes, true, wbSol, Var, Data, Config, DistrInfo, Fixed);

	status = [status_Tras, ' ', status_Dfig, ' ', status_Pv, ' ', status_ClNI, ' ', status_ClRes, ' ', status_Alm, ' ', status_OpSis];

end

function [Var_Prob, opt_Prob, opt_Acum, status_Prob_i] = solve_problem(strProb, previo, totalProb, opt_Prob, opt_Acum, problem, it, NodesSet, oneTime, wbSol, Var_Prev, Data_d, Config, DistrInfo, conFix)

	status_Prob_i = [];
	if ~isempty(NodesSet)
		DistrInfo.Bus = NodesSet;
		if oneTime
			[Var_Prob, opt_Prob, opt_Acum, status_Prob_i] = solve_one(strProb, previo, totalProb, opt_Prob, opt_Acum, problem, it, NodesSet, wbSol, Var_Prev, Data_d, Config, DistrInfo, conFix);
		else
			Var_Prob = Var_Prev;
			for i = 1:length(NodesSet)
				DistrInfo.Bus = NodesSet(i);
				[Var_Prob, opt_Prob, opt_Acum, status_P] = solve_one(strProb, previo+i-1, totalProb, opt_Prob, opt_Acum, problem, it, NodesSet, wbSol, Var_Prob, Data_d, Config, DistrInfo, conFix);
				status_Prob_i = [status_Prob_i num2str(i) ': ' status_P ' '];
			end
		end
	else
		Var_Prob = Var_Prev;
	end
	leyenda = [strProb ' - Finished' ' - ' status_Prob_i]
end

function [Var_Prob, opt_Prob, opt_Acum, status_Prob_i] = solve_one(strProb, previo, totalProb, opt_Prob, opt_Acum, problem, it, NodesSet, wbSol, Var_Prev, Data_d, Config, DistrInfo, conFix)

% 	waitbar(previo/(totalProb),wbSol,sprintf(['Calculo de solucion distr - it ', num2str(it), ' - ', strProb]))
		[Var_Prob_i, opt_Aux, status_Prob_i] = problem(Data_d, Config, DistrInfo);
% 	waitbar(previo+1/(totalProb),wbSol,sprintf(['Calculo de solucion distr - it ', num2str(it), ' - ', strProb]))

	Var_Prob = fuseStructs(Var_Prev, Var_Prob_i);
	opt_Acum = opt_Acum + opt_Aux;
	opt_Prob(it,:) = opt_Prob(it,:) + opt_Aux;
end

function [ret] = isNotSolved(status)
 ret = strcmp(status, 'Unbounded') || strcmp(status, 'Infeasible') || strcmp(status, 'Failed') || strcmp(status, 'Overdetermined');
end
% if ~strcmp(status, 'Solved') && ~strcmp(status, 'Suboptimal')
% 	save(Config.workspace_var_file);
% end