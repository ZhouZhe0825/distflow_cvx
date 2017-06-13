%% Declaracion de parametros de simulacion

DflowD = cargaDflowDDefault();

% Trafos
Trafo1.N = [-2 2];
Trafo1.TP = .005;
Trafo1.nodI = 1;
Trafo1.nodJ = 2;
Trafo1.ini = 0;
Trafo1.cambio = 1;
DflowD.Trafos = [Trafo1];

% Caps
Cap1.N = [0 3];
Cap1.TP = .005;
Cap1.nod = 9;
Cap1.ini = 0;
Cap1.cambio = 1;
DflowD.Caps = [Cap1];

% Cargas
Carga1.pC = 0.5;
Carga1.qC = Carga1.pC *.015;
Carga1.dur = 3;
Carga1.nMultipTop = 1.1;
Carga1.nMultipLow = .75;
Carga1.fileU = 'casos\util\betaE.csv';
Carga1.nod = 5;

Carga2.pC = 0.25;
Carga2.qC = Carga2.pC *.015;
Carga2.dur = 2;
Carga2.nMultipTop = 1.1;
Carga2.nMultipLow = .75;
Carga2.fileU = 'casos\util\betaE.csv';
Carga2.nod = 4;

DflowD.Cargas = [Carga1;Carga2];

%% Nombres de archivos
% 
DflowD.inFilename = 'casos\PU_example\PU_example4.xls';
DflowD.fileCurvaCarga = 'casos\PU_example\cargas\carga_PU_example.csv';
DflowD.fileUtilBetaT = 'casos\util\betaT.csv';
DflowD.fileTemp = 'casos\temp\tempInvierno.csv';
DflowD.fileCostosTension = 'casos\costos\tension\costosTension.csv';
DflowD.fileCostosTras = 'casos\costos\trasmision\costosTrasmision.csv';

%% Configuracion de simulacion
iniEtapa = 1;
CantHorasEtapa = 1;

Config.iniEtapa = iniEtapa;
Config.Etapas = 4*CantHorasEtapa;
Config.outFilename = 'PU_example4_tap_cap_ch';
Config.runNxN = true;
Config.runM = true;

% Centralizado
Config.Centr = [];
% Mosek
% Config.Centr{1, 1} = 'MSK_DPAR_MIO_TOL_REL_RELAX_INT';
% Config.Centr{1, 2} = 0.01;
% Config.Centr{2, 1} = 'MSK_DPAR_OPTIMIZER_MAX_TIME';
% Config.Centr{2, 2} = 600;
% Gurobi
% Config.Centr{1, 1} = 'MIPGap';
% Config.Centr{1, 2} = 0.01;1
% Config.Centr{2, 1} = 'TimeLimit';
% Config.Centr{2, 2} = 600;

%% Carga de datos

[Data] = loadData(DflowD);

%% Correr Simulacion

runSimulation(Data, DflowD, Config);

