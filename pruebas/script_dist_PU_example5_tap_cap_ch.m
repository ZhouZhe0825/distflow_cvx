%% Declaracion de parametros de simulacion

DflowD = cargaDflowDDefault();

% Trafos
DflowD.Trafos(1,1).N = [-8 8];
DflowD.Trafos(1,1).TP = .025;
DflowD.Trafos(1,1).nodI = 1;
DflowD.Trafos(1,1).nodJ = 2;
DflowD.Trafos(1,1).ini = 0;
DflowD.Trafos(1,1).reg = 0;
DflowD.Trafos(1,1).cambio = 1;

% Caps
DflowD.Cap(1,1).N = [0 3];
DflowD.Cap(1,1).TP = .005;
DflowD.Cap(1,1).nod = 9;
DflowD.Cap(1,1).ini = 0;
DflowD.Cap(1,1).cambio = 1;

% Cargas
DflowD.Cargas(1,1).pC = 0.5;
DflowD.Cargas(1,1).qC = DflowD.Cargas(1,1).pC *.015;
DflowD.Cargas(1,1).dur = 4;
DflowD.Cargas(1,1).nMultipTop = 1.1;
DflowD.Cargas(1,1).nMultipLow = .75;
DflowD.Cargas(1,1).fileU = 'casos\util\betaE.csv';
DflowD.Cargas(1,1).nod = 5;

DflowD.Cargas(2,1).pC = 0.25;
DflowD.Cargas(2,1).qC = DflowD.Cargas(2,1).pC *.015;
DflowD.Cargas(2,1).dur = 2;
DflowD.Cargas(2,1).nMultipTop = 1.1;
DflowD.Cargas(2,1).nMultipLow = .75;
DflowD.Cargas(2,1).fileU = 'casos\util\betaE.csv';
DflowD.Cargas(2,1).nod = 9;

% Trasmision
DflowD.Tras.nod = 2;
DflowD.Tras.fileG = 'casos\gen\tras\tras_def.csv';
DflowD.Tras.uLow = 1;
DflowD.Tras.uTop = 1;
DflowD.Tras.fileC = 'casos\costos\trasmision\costosTrasmision.csv';

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
Config.outFilename = 'PU_example5_tap_cap_ch_dist';

Config.Distr.gama_ini = .25;
Config.Distr.gama_fin = .25;
Config.Distr.TopIT = 6;
Config.Distr.tol = 5e-8;


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

Config.SubP{1, 1} = 'MSK_DPAR_OPTIMIZER_MAX_TIME';
Config.SubP{1, 2} = 172800;
Config.SubP{2, 1} = 'MSK_DPAR_MIO_TOL_REL_GAP';
Config.SubP{2, 2} = 1;

Config.Ini{1, 1} = 'MSK_DPAR_OPTIMIZER_MAX_TIME';
Config.Ini{1, 2} = 172800;
Config.Ini{2, 1} = 'MSK_DPAR_MIO_TOL_REL_GAP';
Config.Ini{2, 2} = 1;


%% Carga de datos

[Data] = loadData(DflowD);

%% Correr Simulacion

runSimulation(Data, DflowD, Config);

