function runSimulation(Data, Trafos, Caps, Switches, App, Eolicos, Solares, Baterias, Cargas, ACs, inFilename, fileCurvaCarga, fileUtilBetaT, fileTemp, fileCostosTension, fileCostosTras, Config)
	cantTaps = length(Trafos);
	cantCaps = length(Caps);
	cantCargs = length(Cargas);

	date_post = datestr(now,'yyyymmdd_HHMMSS');

	outputDir = ['resultados\',Config.outFilename,'\',date_post];
	outputDirInputs = [outputDir,'\','entradas'];
	outputDirOutputs = [outputDir,'\','salidas'];
	mkdir(outputDirInputs);
	mkdir(outputDirOutputs);

	for i = 1:length(Eolicos)
		eol = Eolicos(i);
		newFileG = ['eol',num2str(i),'fileG.csv'];
		newFileC = ['eol',num2str(i),'fileC.csv'];
		copyrename(eol.fileG,newFileG,outputDirInputs);
		copyrename(eol.fileC,newFileC,outputDirInputs);
		eol.fileG = [outputDirInputs, '\', newFileG];
		eol.fileC = [outputDirInputs, '\', newFileC];
	end
	for i = 1:length(Solares)
		sol = Solares(i);
		newFileG = ['sol',num2str(i),'fileG.csv'];
		newFileC = ['sol',num2str(i),'fileC.csv'];
		copyrename(sol.fileG,newFileG,outputDirInputs);
		copyrename(sol.fileC,newFileC,outputDirInputs);
		sol.fileG = [outputDirInputs, '\', newFileG];
		sol.fileC = [outputDirInputs, '\', newFileC];
	end
	for i = 1:length(Cargas)
		carg = Cargas(i);
		newFileU = ['carg',num2str(i),'fileU.csv'];
		copyrename(carg.fileU,newFileU,outputDirInputs);
		carg.fileU = [outputDirInputs, '\', newFileU];
	end
	for i = 1:length(ACs)
		ac = ACs(i);
		newFileT = ['ac',num2str(i),'fileT.csv'];
		copyrename(ac.fileT,newFileT,outputDirInputs);
		ac.fileT = [outputDirInputs, '\', newFileT];
	end
	copyrename(inFilename,'inFilename.xls',outputDirInputs);
	inFilename = [outputDirInputs, '\inFilename.xls'];

	copyrename(fileCurvaCarga,'fileCurvaCarga.csv',outputDirInputs);
	fileCurvaCarga = [outputDirInputs, '\fileCurvaCarga.csv'];

	copyrename(fileUtilBetaT,'betaT.csv',outputDirInputs);
	fileUtilBetaT = [outputDirInputs, '\betaT.csv'];

	copyrename(fileTemp,'fileTemp.csv',outputDirInputs);
	fileTemp = [outputDirInputs, '\fileTemp.csv'];

	copyrename(fileCostosTension,'fileCostosTension.csv',outputDirInputs);
	fileCostosTension = [outputDirInputs, '\fileCostosTension.csv'];

	copyrename(fileCostosTras,'fileCostosTras.csv',outputDirInputs);
	fileCostosTras = [outputDirInputs, '\fileCostosTras.csv'];

	save([outputDirInputs,'\inputs.mat'],...
		'Data','inFilename','Trafos','Caps','Cargas','App','Switches','fileCurvaCarga','Eolicos','Solares',...
		'fileUtilBetaT','fileTemp','Baterias', 'ACs','fileCostosTension','fileCostosTras','Config');

	diary([outputDirOutputs, '\console.log']);

	% copyfile('.gitignore', 'alsdjfad')
	% exist('alsdjfad')
	% diary('prueba.log')
	% diary('off')


	leyenda = ['------------------------------------ ' Config.outFilename ' ------------------------------------']

    Var_nxn = [];
    opt_nxn = [];
    DataNxN = [];
    if Config.runNxN
    	[Var_nxn, opt_nxn, DataNxN] = llamarCentralizadoNxN(Data, Config);
    else
        [DataNxN] = reshapeDataNxN(Data, Config);
    end

    Var_m = [];
    opt_m = [];
    DataM = [];
    if Config.runM
    	[Var_m, opt_m, DataM] = llamarCentralizadoM(Data, Config);
    end

	save([outputDirOutputs, '\outputs.mat'],'Var_nxn','opt_nxn','DataNxN','Config','Var_m','opt_m','DataM');

	diary('off');

    if Config.runNxN
        printSalidasDistflowNxN(Var_nxn, DataNxN, Config, cantTaps, cantCaps, cantCargs, [outputDirOutputs, '\output_nxn'], [], [], [], [], []);
    end
    if Config.runM
        printSalidasDistflowM(Var_m, DataNxN, Config, cantTaps, cantCaps, cantCargs, [outputDirOutputs, '\output_m'], [], [], [], [], []);
    end
	
end