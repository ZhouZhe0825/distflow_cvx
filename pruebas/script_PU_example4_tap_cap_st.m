%% Declaracion de parametros de simulacion

DflowD = cargaDflowDDefault();

% Baterias
DflowD.Baterias(1,1).nod = 5;
DflowD.Baterias(1,1).type = 'casos\bat\bat_def.csv';
DflowD.Baterias(1,1).EIni = .5;
DflowD.Baterias(1,1).pgIni = 0;
DflowD.Baterias(1,1).qgIni = 0;

DflowD.Baterias(2,1).nod = 8;
DflowD.Baterias(2,1).type = 'casos\bat\bat_def.csv';
DflowD.Baterias(2,1).EIni = .5;
DflowD.Baterias(2,1).pgIni = 0;
DflowD.Baterias(2,1).qgIni = 0;
DflowD.Baterias(2,1).pgIni = 0;
DflowD.Baterias(2,1).qgIni = 0;

% Trafos
DflowD.Trafos(1,1).N = [-8 8];
DflowD.Trafos(1,1).TP = .025;
DflowD.Trafos(1,1).nodI = 1;
DflowD.Trafos(1,1).nodJ = 2;
DflowD.Trafos(1,1).ini = 0;
DflowD.Trafos(1,1).reg = 0;
DflowD.Trafos(1,1).cambio = 1;

% Caps
DflowD.Caps(1,1).N = [0 3];
DflowD.Caps(1,1).TP = .005;
DflowD.Caps(1,1).nod = 9;
DflowD.Caps(1,1).ini = 0;
DflowD.Caps(1,1).cambio = 1;

% Trasmision
DflowD.Tras(1,1).nod = 1;
DflowD.Tras(1,1).fileG = 'casos\gen\tras\tras_def.csv';
DflowD.Tras(1,1).uLow = 1;
DflowD.Tras(1,1).uTop = 1;
DflowD.Tras(1,1).fileC = 'casos\costos\trasmision\costosTrasmision.csv';
DflowD.Tras(1,1).pgIni = 0;
DflowD.Tras(1,1).qgIni = 0;

%% Nombres de archivos
% 
DflowD.inFilename = 'casos\PU_example\PU_example4.xls';
DflowD.fileCurvaCarga = 'casos\PU_example\cargas\carga_PU_example.csv';
DflowD.fileUtilBetaT = 'casos\util\betaT.csv';
DflowD.utilOptFuncCuad = true;
DflowD.fileTemp = 'casos\temp\tempInvierno.csv';
DflowD.fileCostosTension = 'casos\costos\tension\costosTension.csv';

%% Configuracion de simulacion
iniEtapa = 1;
CantHorasEtapa = 1;

Config = [];
Config.iniEtapa = iniEtapa;
Config.Etapas = 4*CantHorasEtapa;
Config.outFilename = 'PU_example4_tap_cap_st';

% Centralizado
Config.Centr = [];
% Mosek
Config.Centr{1, 1} = 'MSK_DPAR_MIO_TOL_REL_RELAX_INT';
Config.Centr{1, 2} = 0.01;
Config.Centr{2, 1} = 'MSK_DPAR_OPTIMIZER_MAX_TIME';
Config.Centr{2, 2} = 120;
% Gurobi
% Config.Centr{1, 1} = 'MIPGap';
% Config.Centr{1, 2} = 0.01;1
% Config.Centr{2, 1} = 'TimeLimit';
% Config.Centr{2, 2} = 600;

%% Carga de datos

[Data] = loadData(DflowD);

%% Correr Simulacion

runSimulation(Data, DflowD, Config);

