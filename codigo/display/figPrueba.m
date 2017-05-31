function figPrueba()


mainW.fig = figure('Name', 'Distflow - Ort', 'MenuBar', 'none', 'NumberTitle', 'off');
mainW.pan_datosSim = uipanel(mainW.fig,'Title','Datos Simulacion','Position',[0 .1 .5 .8]);
mainW.pan_redSim = uipanel(mainW.fig,'Title','Red a Simular','Position',[.5 .1 .5 .8]); 
mainW.tabgr_Carga = uitabgroup('Parent',mainW.pan_datosSim,'Position',[0 .15 1 .8]);
mainW.tab_General = uitab('Parent',mainW.tabgr_Carga,'Title','General');
mainW.axes_dispRed = axes('Parent',mainW.pan_redSim,'Position',[.05,.05,.9,.9], 'XTick', [], 'YTick', []);
mainW.but_carga = uicontrol(mainW.pan_datosSim,'Style','pushbutton', 'String', 'Carga', 'Units', 'normalized', 'Position', [.025 .025 .25 .1],'Callback',@but_carga_callback);
mainW.but_ejecutar = uicontrol(mainW.pan_datosSim,'Style','pushbutton', 'String', 'Ejecutar', 'Units', 'normalized', 'Position', [.725 .025 .25 .1],'Callback',@but_ejecutar_callback);

mainW.red_file = uicontrol(mainW.tab_General,'Style','edit', 'Units', 'normalized', 'Position', [.5 .65 .4 .1]);
mainW.label_red_file = uicontrol(mainW.tab_General,'Style','text', 'String', 'Archivo de red', 'Units', 'normalized', 'Position', [.1 .65 .4 .1]);

mainW.carga_file = uicontrol(mainW.tab_General,'Style','edit', 'Units', 'normalized', 'Position', [.5 .5 .4 .1]);
mainW.label_carga_file = uicontrol(mainW.tab_General,'Style','text', 'String', 'Archivo de carga', 'Units', 'normalized', 'Position', [.1 .5 .4 .1]);

mainW.betaT_file = uicontrol(mainW.tab_General,'Style','edit', 'Units', 'normalized', 'Position', [.5 .35 .4 .1]);
mainW.label_betaT_file = uicontrol(mainW.tab_General,'Style','text', 'String', 'Archivo de betaT', 'Units', 'normalized', 'Position', [.1 .35 .4 .1]);

mainW.temp_file = uicontrol(mainW.tab_General,'Style','edit', 'Units', 'normalized', 'Position', [.5 .2 .4 .1]);
mainW.label_temp_file = uicontrol(mainW.tab_General,'Style','text', 'String', 'Archivo de temperatura', 'Units', 'normalized', 'Position', [.1 .2 .4 .1]);

mainW.costosTension_file = uicontrol(mainW.tab_General,'Style','edit', 'Units', 'normalized', 'Position', [.5 .05 .4 .1]);
mainW.label_costosTension_file = uicontrol(mainW.tab_General,'Style','text', 'String', 'Archivo de costosV', 'Units', 'normalized', 'Position', [.1 .05 .4 .1]);


mainW.Context.Data = [];
mainW.Context.Trafos = [];
mainW.Context.ACs = [];

mainW.Context.App(1).I = 1;
mainW.Context.App(1).Pref = .8;
mainW.Context.App(1).Low = .8;
mainW.Context.App(1).Top = .8;
mainW.Context.App(1).nMultipTop = 1.05;
mainW.Context.App(1).nMultipLow = .95;
mainW.Context.App(1).alpha = 0;
mainW.Context.App(1).tgPhi = .2;
mainW.Context.App(2).I = 2;
mainW.Context.App(2).Pref = .2;
mainW.Context.App(2).Low = 0;
mainW.Context.App(2).Top = .2;
mainW.Context.App(2).nMultipTop = 1.05;
mainW.Context.App(2).nMultipLow = .95;
mainW.Context.App(2).alpha = 0.1;
mainW.Context.App(2).tgPhi = .2;

mainW.Context.Baterias = [];
mainW.Context.Caps = [];
mainW.Context.Cargas = [];
mainW.Context.Eolicos = [];
mainW.Context.fileCostosTension = [];
mainW.Context.fileCostosTras = [];
mainW.Context.fileCurvaCarga = [];
mainW.Context.fileTemp = [];
mainW.Context.fileUtilBetaT = [];
mainW.Context.inFilename = [];
mainW.Context.Solares = [];

mainW.Context.Switches.all = true;
mainW.Context.Switches.i = [];
mainW.Context.Switches.j = [];
mainW.Context.Switches.cY = .1;


mainW.Context.Config.iniEtapa = 1;
mainW.Context.Config.Etapas = 4;
mainW.Context.Config.outFilename = 'PU_example4';
mainW.Context.Config.runNxN = true;
mainW.Context.Config.runM = true;
mainW.Context.Config.Centr = [];


    function but_carga_callback(source,eventdata)
		
		mainW.Context.inFilename = get(mainW.red_file, 'string');
		mainW.Context.fileCurvaCarga = get(mainW.carga_file, 'string');
		mainW.Context.fileUtilBetaT = get(mainW.betaT_file, 'string');
		mainW.Context.fileTemp = get(mainW.temp_file, 'string');
		mainW.Context.fileCostosTension = get(mainW.costosTension_file, 'string');
		mainW.Context.fileCostosTras = 'casos\costos\trasmision\costosTrasmision.csv';

% 		mainW.Context.inFilename = 'casos\PU_example\PU_example4.xls';
% 		mainW.Context.fileCurvaCarga = 'casos\PU_example\cargas\carga_PU_example.csv';
% 		mainW.Context.fileUtilBetaT = 'casos\util\betaT.csv';
% 		mainW.Context.fileTemp = 'casos\temp\tempInvierno.csv';
% 		mainW.Context.fileCostosTension = 'casos\costos\tension\costosTension.csv';
% 		mainW.Context.fileCostosTras = 'casos\costos\trasmision\costosTrasmision.csv';
		

		[mainW.Context.Data] = loadData(...
			mainW.Context.Trafos, mainW.Context.Caps, mainW.Context.Switches, ...
			mainW.Context.App, mainW.Context.Eolicos, mainW.Context.Solares, ...
			mainW.Context.Baterias, mainW.Context.Cargas, mainW.Context.ACs, ...
			mainW.Context.inFilename, mainW.Context.fileCurvaCarga, mainW.Context.fileUtilBetaT, ...
			mainW.Context.fileTemp, mainW.Context.fileCostosTension, mainW.Context.fileCostosTras ...
			);
		plotRed(mainW.Context.Data);
        
    end
	
	function but_ejecutar_callback(source,eventdata)
	
		runSimulation(...
			mainW.Context.Data, mainW.Context.Trafos, mainW.Context.Caps, mainW.Context.Switches, ...
			mainW.Context.App, mainW.Context.Eolicos, mainW.Context.Solares, mainW.Context.Baterias, ...
			mainW.Context.Cargas, mainW.Context.ACs, mainW.Context.inFilename, mainW.Context.fileCurvaCarga, ...
			mainW.Context.fileUtilBetaT, mainW.Context.fileTemp, mainW.Context.fileCostosTension, ...
			mainW.Context.fileCostosTras, mainW.Context.Config ...
			);

	end
end