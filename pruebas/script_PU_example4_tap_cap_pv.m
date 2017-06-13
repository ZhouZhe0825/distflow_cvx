%% Declaracion de parametros de simulacion

DflowD = cargaDflowDDefault();

% Solar
Pv1.nod = 4;
Pv1.type = 'casos\gen\pv\pv_def.csv';
Pv1.fileG = 'casos\gen\pv\pvgen.csv';
Pv1.fileC = 'casos\costos\pv\costosPv.csv';

Pv2.nod = 7;
Pv2.type = 'casos\gen\pv\pv_def.csv';
Pv2.fileG = 'casos\gen\pv\pvgen.csv';
Pv2.fileC = 'casos\costos\pv\costosPv.csv';

DflowD.Solares = [Pv1; Pv2];

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
Config.outFilename = 'PU_example4_tap_cap_pv';
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

