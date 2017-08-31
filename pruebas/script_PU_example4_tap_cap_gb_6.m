%% Declaracion de parametros de simulacion

DflowD = cargaDflowDDefault();

DflowD.App(1,1).I = 1;
DflowD.App(1,1).Pref = 1;
DflowD.App(1,1).Low = 1;
DflowD.App(1,1).Top = 1;
DflowD.App(1,1).nMultipTop = 1.05;
DflowD.App(1,1).nMultipLow = .95;
DflowD.App(1,1).alpha = 0.4;
DflowD.App(1,1).tgPhi = .2;
DflowD.App(2,1).I = 2;
DflowD.App(2,1).Pref = 0;
DflowD.App(2,1).Low = 0;
DflowD.App(2,1).Top = 0;
DflowD.App(2,1).nMultipTop = 0;
DflowD.App(2,1).nMultipLow = 0;
DflowD.App(2,1).alpha = 0;
DflowD.App(2,1).tgPhi = 0;

% Generador Basico
DflowD.GenBas(1,1).nod = 7;
DflowD.GenBas(1,1).fileG = 'casos\gen\basic\genbas_fix.csv';
DflowD.GenBas(1,1).fileC = 'casos\costos\basic\costosBasic_fix.csv';

% Trafos
DflowD.Trafos(1,1).N = [-8 8];
DflowD.Trafos(1,1).TP = .025;
DflowD.Trafos(1,1).nodI = 1;
DflowD.Trafos(1,1).nodJ = 2;
DflowD.Trafos(1,1).ini = 0;
DflowD.Trafos(1,1).reg = 0;
DflowD.Trafos(1,1).cambio = 1;

% Caps
% Trasmision
DflowD.Tras.nod = 1;
DflowD.Tras.fileG = 'casos\gen\tras\tras_def.csv';
DflowD.Tras.uLow = 1;
DflowD.Tras.uTop = 1;
DflowD.Tras.fileC = 'casos\costos\trasmision\costosTrasmision_incr5.csv';

% Cargas

% Switches
DflowD.Switches.i = [5 7];
DflowD.Switches.j = [6 8];
DflowD.Switches.cY = 10;
DflowD.Switches.all = false;

% Baterias
DflowD.Baterias(1,1).nod = 5;
DflowD.Baterias(1,1).type = 'casos\bat\bat_eta1.csv';
DflowD.Baterias(1,1).EIni = .1;

%% Nombres de archivos
% 
DflowD.inFilename = 'casos\PU_example\PU_example4_uTL_aj.xls';
DflowD.fileCurvaCarga = 'casos\PU_example\cargas\carga_PU_example_variable.csv';
DflowD.fileUtilBetaT = 'casos\util\betaT.csv';
DflowD.utilOptFuncCuad = false;
DflowD.fileTemp = 'casos\temp\tempInvierno.csv';
DflowD.fileCostosTension = 'casos\costos\tension\costosTension.csv';

%% Configuracion de simulacion
iniEtapa = 1;
CantHorasEtapa = 24;

Config = [];
Config.iniEtapa = iniEtapa;
Config.Etapas = 4*CantHorasEtapa;
Config.outFilename = 'PU_example4_tap_cap_gb_6';
Config.runNxN = false;
Config.runM = true;

% Centralizado
Config.Centr = [];
% Mosek
Config.Centr{1, 1} = 'MSK_DPAR_MIO_NEAR_TOL_REL_GAP';
Config.Centr{1, 2} = 0.5;
Config.Centr{2, 1} = 'MSK_DPAR_OPTIMIZER_MAX_TIME';
Config.Centr{2, 2} = 300;
% Gurobi
% Config.Centr{1, 1} = 'MIPGap';
% Config.Centr{1, 2} = 0.01;1
% Config.Centr{2, 1} = 'TimeLimit';
% Config.Centr{2, 2} = 600;

%% Carga de datos

[Data] = loadData(DflowD);

%% Correr Simulacion

runSimulation(Data, DflowD, Config);

