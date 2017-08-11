function DflowD = cargaDflowDDefault()

	DflowD.Trafos = [];
	DflowD.ACs = [];

	DflowD.App(1,1).I = 1;
	DflowD.App(1,1).Pref = .8;
	DflowD.App(1,1).Low = .8;
	DflowD.App(1,1).Top = .8;
	DflowD.App(1,1).nMultipTop = 1.05;
	DflowD.App(1,1).nMultipLow = .95;
	DflowD.App(1,1).alpha = 0;
	DflowD.App(1,1).tgPhi = .2;
	DflowD.App(2,1).I = 2;
	DflowD.App(2,1).Pref = .2;
	DflowD.App(2,1).Low = 0;
	DflowD.App(2,1).Top = .2;
	DflowD.App(2,1).nMultipTop = 1.05;
	DflowD.App(2,1).nMultipLow = .95;
	DflowD.App(2,1).alpha = 0.1;
	DflowD.App(2,1).tgPhi = .2;

	DflowD.Baterias = [];
	DflowD.Caps = [];
	DflowD.Cargas = [];
	DflowD.Eolicos = [];
    DflowD.GenBas = [];
	DflowD.fileCostosTension = '';
	DflowD.fileCostosTras = '';
	DflowD.fileCurvaCarga = '';
	DflowD.fileTemp = '';
	DflowD.fileUtilBetaT = '';
	DflowD.inFilename = '';
	DflowD.Solares = [];

	DflowD.Switches.all = true;
	DflowD.Switches.i = [];
	DflowD.Switches.j = [];
	DflowD.Switches.cY = 0;

end