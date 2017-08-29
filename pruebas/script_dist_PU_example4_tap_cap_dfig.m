%% Declaracion de parametros de simulacion

DflowD = cargaDflowDDefault();

% Eolicos
DflowD.Eolicos(1,1).nod = 5;
DflowD.Eolicos(1,1).type = 'casos\gen\dfig\Dfig_1mw_def.csv';
DflowD.Eolicos(1,1).fileG = 'casos\gen\dfig\Dfig_1mw_P_n_.csv';
DflowD.Eolicos(1,1).fileC = 'casos\costos\dfig\costosDfig.csv';


% Trafos
DflowD.Trafos(1,1).N = [-8 8];
DflowD.Trafos(1,1).TP = .025;
DflowD.Trafos(1,1).nodI = 1;
DflowD.Trafos(1,1).nodJ = 2;
DflowD.Trafos(1,1).ini = 0;
DflowD.Trafos(1,1).reg = 0;
DflowD.Trafos(1,1).cambio = 1;

DflowD.Trafos(2,1).N = [-1 1];
DflowD.Trafos(2,1).TP = .1;
DflowD.Trafos(2,1).nodI = 4;
DflowD.Trafos(2,1).nodJ = 5;
DflowD.Trafos(2,1).ini = 0;
DflowD.Trafos(2,1).reg = 1;
DflowD.Trafos(2,1).cambio = .1;

% Caps
DflowD.Cap(1,1).N = [0 3];
DflowD.Cap(1,1).TP = .005;
DflowD.Cap(1,1).nod = 9;
DflowD.Cap(1,1).ini = 0;
DflowD.Cap(1,1).cambio = 1;

%% Nombres de archivos
% 
DflowD.inFilename = 'casos\PU_example\PU_example4.xls';
DflowD.fileCurvaCarga = 'casos\PU_example\cargas\carga_PU_example.csv';
DflowD.fileUtilBetaT = 'casos\util\betaT.csv';
DflowD.utilOptFuncCuad = true;
DflowD.fileTemp = 'casos\temp\tempInvierno.csv';
DflowD.fileCostosTension = 'casos\costos\tension\costosTension.csv';
DflowD.fileCostosTras = 'casos\costos\trasmision\costosTrasmision.csv';

%% Configuracion de simulacion
iniEtapa = 1;
CantHorasEtapa = 1;

Config.iniEtapa = iniEtapa;
Config.Etapas = 4*CantHorasEtapa;
Config.outFilename = 'PU_example4_tap_cap_dfig_dist';
Config.runNxN = false;
Config.runM = true;

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

