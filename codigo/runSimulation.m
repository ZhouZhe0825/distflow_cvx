function runSimulation(Data, DflowD, Config)
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
	for i = 1:length(DflowD.GenBas)
		gbas = DflowD.GenBas(i);
		newFileG = ['gbas',num2str(i),'fileG.csv'];
		newFileC = ['gbas',num2str(i),'fileC.csv'];
		copyrename(gbas.fileG,newFileG,outputDirInputs);
		copyrename(gbas.fileC,newFileC,outputDirInputs);
		gbas.fileG = [outputDirInputs, '\', newFileG];
		gbas.fileC = [outputDirInputs, '\', newFileC];
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
	for i = 1:length(DflowD.Tras)
		tra = DflowD.Tras(i);
		newFileG = ['tra',num2str(i),'fileG.csv'];
		newFileC = ['tra',num2str(i),'fileC.csv'];
		copyrename(tra.fileG,newFileG,outputDirInputs);
		copyrename(tra.fileC,newFileC,outputDirInputs);
		tra.fileG = [outputDirInputs, '\', newFileG];
		tra.fileC = [outputDirInputs, '\', newFileC];
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

	save([outputDirInputs,'\inputs.mat'],...
		'Data','DflowD','Config');

	diary([outputDirOutputs, '\console.log']);

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
		save([outputDirOutputs, '\outputs.mat'],'Config','Var_m','opt_m','DataM');

		diary('off');
	else
		[Var_dist_conE, Var_centr, Var_F, opt_dist_conE, opt_centr, opt_F, status, DataM, Ev] = llamarDistribuidoM(Data, Config, Var_centr, opt_centr, Var_ini);

		save([outputDirOutputs, '\outputs.mat'],'Var_dist_conE', 'Var_centr', 'Var_F', 'Var_ini', 'opt_dist_conE', 'opt_centr', 'opt_F', 'Ev','Config','Var_m','opt_m','DataM');

		diary('off');
	end

	if ~isfield(Config,'Distr')
		printSalidasDistflowM(Var_m, DataM, Config, [outputDirOutputs, '\output_m2'], [], [], [], [], []);
	else
        printSalidasDistflowM(Var_F,         DataM, Config, [outputDirOutputs, '\output_m'],           Ev.opt, Ev.mu, Ev.lambda, Ev.difP, Ev.difQ);
		printSalidasDistflowM(Var_centr,     DataM, Config, [outputDirOutputs, '\output_centr'],       [],    [],   [],       [],     []);
		printSalidasDistflowM(Var_dist_conE, DataM, Config, [outputDirOutputs, '\output_dist_conE'],   [],    [],   [],       [],     []);
	end
	
end