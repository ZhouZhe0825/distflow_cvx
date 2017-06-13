function runSimulation(Data, DflowD, Config)
	cantTaps = length(DflowD.Trafos);
	cantCaps = length(DflowD.Caps);
	cantCargs = length(DflowD.Cargas);

	date_post = datestr(now,'yyyymmdd_HHMMSS');

	outputDir = ['resultados\',Config.outFilename,'\',date_post];
	outputDirInputs = [outputDir,'\','entradas'];
	outputDirOutputs = [outputDir,'\','salidas'];
	mkdir(outputDirInputs);
	mkdir(outputDirOutputs);

	for i = 1:length(DflowD.Eolicos)
		eol = DflowD.Eolicos(i);
		newFileG = ['eol',num2str(i),'fileG.csv'];
		newFileC = ['eol',num2str(i),'fileC.csv'];
		copyrename(eol.fileG,newFileG,outputDirInputs);
		copyrename(eol.fileC,newFileC,outputDirInputs);
		eol.fileG = [outputDirInputs, '\', newFileG];
		eol.fileC = [outputDirInputs, '\', newFileC];
	end
	for i = 1:length(DflowD.Solares)
		sol = DflowD.Solares(i);
		newFileG = ['sol',num2str(i),'fileG.csv'];
		newFileC = ['sol',num2str(i),'fileC.csv'];
		copyrename(sol.fileG,newFileG,outputDirInputs);
		copyrename(sol.fileC,newFileC,outputDirInputs);
		sol.fileG = [outputDirInputs, '\', newFileG];
		sol.fileC = [outputDirInputs, '\', newFileC];
	end
	for i = 1:length(DflowD.Cargas)
		carg = DflowD.Cargas(i);
		newFileU = ['carg',num2str(i),'fileU.csv'];
		copyrename(carg.fileU,newFileU,outputDirInputs);
		carg.fileU = [outputDirInputs, '\', newFileU];
	end
	for i = 1:length(DflowD.ACs)
		ac = DflowD.ACs(i);
		newFileT = ['ac',num2str(i),'fileT.csv'];
		copyrename(ac.fileT,newFileT,outputDirInputs);
		ac.fileT = [outputDirInputs, '\', newFileT];
	end
	copyrename(DflowD.inFilename,'inFilename.xls',outputDirInputs);
	DflowD.inFilename = [outputDirInputs, '\inFilename.xls'];

	copyrename(DflowD.fileCurvaCarga,'fileCurvaCarga.csv',outputDirInputs);
	DflowD.fileCurvaCarga = [outputDirInputs, '\fileCurvaCarga.csv'];

	copyrename(DflowD.fileUtilBetaT,'betaT.csv',outputDirInputs);
	DflowD.fileUtilBetaT = [outputDirInputs, '\betaT.csv'];

	copyrename(DflowD.fileTemp,'fileTemp.csv',outputDirInputs);
	DflowD.fileTemp = [outputDirInputs, '\fileTemp.csv'];

	copyrename(DflowD.fileCostosTension,'fileCostosTension.csv',outputDirInputs);
	DflowD.fileCostosTension = [outputDirInputs, '\fileCostosTension.csv'];

	copyrename(DflowD.fileCostosTras,'fileCostosTras.csv',outputDirInputs);
	DflowD.fileCostosTras = [outputDirInputs, '\fileCostosTras.csv'];

	save([outputDirInputs,'\inputs.mat'],...
		'Data','DflowD','Config');

	diary([outputDirOutputs, '\console.log']);

	% copyfile('.gitignore', 'alsdjfad')
	% exist('alsdjfad')
	% diary('prueba.log')
	% diary('off')


	leyenda = ['------------------------------------ ' Config.outFilename ' ------------------------------------']

    Var_nxn = [];
    opt_nxn = [];
    DataNxN = [];

	Var_m = [];
	opt_m = [];
	DataM = [];
	
	Var_centr = [];
	opt_centr = [];
	Var_ini = [];
	opt_ini = [];

	
	if ~isfield(Config,'Distr')
		if Config.runNxN
			[Var_nxn, opt_nxn, DataNxN] = llamarCentralizadoNxN(Data, Config);
		else
			[DataNxN] = reshapeDataNxN(Data, Config);
		end

		if Config.runM
			[Var_m, opt_m, DataM] = llamarCentralizadoM(Data, Config);
		end

		save([outputDirOutputs, '\outputs.mat'],'Var_nxn','opt_nxn','DataNxN','Config','Var_m','opt_m','DataM');

		diary('off');
	else
		[DataNxN] = reshapeDataNxN(Data, Config);

		[Var_dist_conE, Var_centr, Var_F, opt_dist_conE, opt_centr, opt_F, status, DataM, Ev] = llamarDistribuidoM(Data, Config, Var_centr, opt_centr, Var_ini);

		save([outputDirOutputs, '\outputs.mat'],'Var_dist_conE', 'Var_centr', 'Var_F', 'Var_ini', 'opt_dist_conE', 'opt_centr', 'opt_F', 'Ev','DataNxN','Config','Var_m','opt_m','DataM');

		diary('off');
	end

	if ~isfield(Config,'Distr')
		if Config.runNxN
			printSalidasDistflowNxN(Var_nxn, DataNxN, Config, cantTaps, cantCaps, cantCargs, [outputDirOutputs, '\output_nxn'], [], [], [], [], []);
		end
		if Config.runM
			printSalidasDistflowM(Var_m, DataNxN, Config, cantTaps, cantCaps, cantCargs, [outputDirOutputs, '\output_m'], [], [], [], [], []);
		end
	else
		printSalidasDistflowM(Var_F,         DataNxN, Config, cantTaps, cantCaps, cantCargs, [outputDirOutputs, '\output_m'],           Ev.opt, Ev.mu, Ev.lambda, Ev.difP, Ev.difQ);
		printSalidasDistflowM(Var_centr,     DataNxN, Config, cantTaps, cantCaps, cantCargs, [outputDirOutputs, '\output_centr'],       [],    [],   [],       [],     []);
		printSalidasDistflowM(Var_dist_conE, DataNxN, Config, cantTaps, cantCaps, cantCargs, [outputDirOutputs, '\output_dist_conE'],   [],    [],   [],       [],     []);
	end
	
end