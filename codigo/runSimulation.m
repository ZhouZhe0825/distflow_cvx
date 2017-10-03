function runSimulation(Data, DflowD, Config)
	date_post = datestr(now,'yyyymmdd_HHMMSS');

	outputDir = ['resultados\',Config.outFilename,'\',date_post];
	Config.outputDirInputs = [outputDir,'\','entradas'];
	Config.outputDirOutputs = [outputDir,'\','salidas'];
	Config.outputDirErrors = [outputDir,'\','errors'];
	mkdir(Config.outputDirInputs);
	mkdir(Config.outputDirOutputs);

	for i = 1:length(DflowD.Eolicos)
		eol = DflowD.Eolicos(i);
		newFileG = ['eol',num2str(i),'fileG.csv'];
		newFileC = ['eol',num2str(i),'fileC.csv'];
		copyrename(eol.fileG,newFileG,Config.outputDirInputs);
		copyrename(eol.fileC,newFileC,Config.outputDirInputs);
		eol.fileG = [Config.outputDirInputs, '\', newFileG];
		eol.fileC = [Config.outputDirInputs, '\', newFileC];
	end
	for i = 1:length(DflowD.GenBas)
		gbas = DflowD.GenBas(i);
		newFileG = ['gbas',num2str(i),'fileG.csv'];
		newFileC = ['gbas',num2str(i),'fileC.csv'];
		copyrename(gbas.fileG,newFileG,Config.outputDirInputs);
		copyrename(gbas.fileC,newFileC,Config.outputDirInputs);
		gbas.fileG = [Config.outputDirInputs, '\', newFileG];
		gbas.fileC = [Config.outputDirInputs, '\', newFileC];
	end
	for i = 1:length(DflowD.Solares)
		sol = DflowD.Solares(i);
		newFileG = ['sol',num2str(i),'fileG.csv'];
		newFileC = ['sol',num2str(i),'fileC.csv'];
		copyrename(sol.fileG,newFileG,Config.outputDirInputs);
		copyrename(sol.fileC,newFileC,Config.outputDirInputs);
		sol.fileG = [Config.outputDirInputs, '\', newFileG];
		sol.fileC = [Config.outputDirInputs, '\', newFileC];
	end
	for i = 1:length(DflowD.Cargas)
		carg = DflowD.Cargas(i);
		newFileU = ['carg',num2str(i),'fileU.csv'];
		copyrename(carg.fileU,newFileU,Config.outputDirInputs);
		carg.fileU = [Config.outputDirInputs, '\', newFileU];
	end
	for i = 1:length(DflowD.ACs)
		ac = DflowD.ACs(i);
		newFileT = ['ac',num2str(i),'fileT.csv'];
		copyrename(ac.fileT,newFileT,Config.outputDirInputs);
		ac.fileT = [Config.outputDirInputs, '\', newFileT];
	end
	for i = 1:length(DflowD.Tras)
		tra = DflowD.Tras(i);
		newFileG = ['tra',num2str(i),'fileG.csv'];
		newFileC = ['tra',num2str(i),'fileC.csv'];
		copyrename(tra.fileG,newFileG,Config.outputDirInputs);
		copyrename(tra.fileC,newFileC,Config.outputDirInputs);
		tra.fileG = [Config.outputDirInputs, '\', newFileG];
		tra.fileC = [Config.outputDirInputs, '\', newFileC];
	end
	copyrename(DflowD.inFilename,'inFilename.xls',Config.outputDirInputs);
	DflowD.inFilename = [Config.outputDirInputs, '\inFilename.xls'];

	copyrename(DflowD.fileCurvaCarga,'fileCurvaCarga.csv',Config.outputDirInputs);
	DflowD.fileCurvaCarga = [Config.outputDirInputs, '\fileCurvaCarga.csv'];

	copyrename(DflowD.fileUtilBetaT,'betaT.csv',Config.outputDirInputs);
	DflowD.fileUtilBetaT = [Config.outputDirInputs, '\betaT.csv'];

	copyrename(DflowD.fileTemp,'fileTemp.csv',Config.outputDirInputs);
	DflowD.fileTemp = [Config.outputDirInputs, '\fileTemp.csv'];

	copyrename(DflowD.fileCostosTension,'fileCostosTension.csv',Config.outputDirInputs);
	DflowD.fileCostosTension = [Config.outputDirInputs, '\fileCostosTension.csv'];

	save([Config.outputDirInputs,'\inputs.mat'],...
		'Data','DflowD','Config');

	diary([Config.outputDirOutputs, '\console.log']);

	% copyfile('.gitignore', 'alsdjfad')
	% exist('alsdjfad')
	% diary('prueba.log')
	% diary('off')


	leyenda = ['------------------------------------ ' Config.outFilename ' ------------------------------------']

	Var_m = [];
	opt_m = [];
	DataM = [];
	
	Var_centr = [];
	opt_centr = [];
	Var_ini = [];
	opt_ini = [];

	
	if ~isfield(Config,'Distr')
		[Var_m, opt_m, DataM] = llamarCentralizadoM(Data, Config);
		save([Config.outputDirOutputs, '\outputs.mat'],'Config','Var_m','opt_m','DataM');

		diary('off');
	else
		[Var_dist_conE, Var_centr, Var_F, Var_ini, opt_dist_conE, opt_centr, opt_F, opt_ini, status, DataM, Ev, error] = llamarDistribuidoM(Data, Config, Var_centr, opt_centr, Var_ini);
		save([Config.outputDirOutputs, '\outputs.mat'],'Var_dist_conE', 'Var_centr', 'Var_F', 'Var_ini', 'opt_dist_conE', 'opt_centr', 'opt_F', 'opt_ini', 'Ev','Config','Var_m','opt_m','DataM');

		diary('off');
	end

	if ~isfield(Config,'Distr')
		printSalidasDistflowM(Var_m, DataM, Config, [Config.outputDirOutputs, '\output_m'], [], [], [], [], []);
	else
		printSalidasDistflowM(Var_F,         DataM, Config, [Config.outputDirOutputs, '\output_m'],           Ev.opt, Ev.mu, Ev.lambda, Ev.difPAbs, Ev.difQAbs);
		printSalidasDistflowM(Var_centr,     DataM, Config, [Config.outputDirOutputs, '\output_centr'],       [],    [],   [],       [],     []);
		printSalidasDistflowM(Var_ini,       DataM, Config, [Config.outputDirOutputs, '\output_ini'],         [],    [],   [],       [],     []);
		printSalidasDistflowM(Var_dist_conE, DataM, Config, [Config.outputDirOutputs, '\output_dist_conE'],   [],    [],   [],       [],     []);
	end
	
end