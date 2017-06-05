function DflowD = cargaDflowDDefault()

	DflowD.Data = [];
	DflowD.Trafos = [];
	DflowD.ACs = [];

	DflowD.App(1).I = 1;
	DflowD.App(1).Pref = .8;
	DflowD.App(1).Low = .8;
	DflowD.App(1).Top = .8;
	DflowD.App(1).nMultipTop = 1.05;
	DflowD.App(1).nMultipLow = .95;
	DflowD.App(1).alpha = 0;
	DflowD.App(1).tgPhi = .2;
	DflowD.App(2).I = 2;
	DflowD.App(2).Pref = .2;
	DflowD.App(2).Low = 0;
	DflowD.App(2).Top = .2;
	DflowD.App(2).nMultipTop = 1.05;
	DflowD.App(2).nMultipLow = .95;
	DflowD.App(2).alpha = 0.1;
	DflowD.App(2).tgPhi = .2;

	DflowD.Baterias = [];
	DflowD.Caps = [];
	DflowD.Cargas = [];
	DflowD.Eolicos = [];
	DflowD.fileCostosTension = [];
	DflowD.fileCostosTras = [];
	DflowD.fileCurvaCarga = [];
	DflowD.fileTemp = [];
	DflowD.fileUtilBetaT = [];
	DflowD.inFilename = [];
	DflowD.Solares = [];

	DflowD.Switches.all = true;
	DflowD.Switches.i = [];
	DflowD.Switches.j = [];
	DflowD.Switches.cY = 0;

end