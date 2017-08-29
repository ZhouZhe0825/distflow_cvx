%% Declaracion de parametros de simulacion

DflowD = cargaDflowDDefault();

% Eolicos
DflowD.Eolicos(1,1).nod = 5;
DflowD.Eolicos(1,1).type = 'casos\gen\dfig\Dfig_1mw_def.csv';
DflowD.Eolicos(1,1).fileG = 'casos\gen\dfig\Dfig_1mw_P_n_.csv';
DflowD.Eolicos(1,1).fileC = 'casos\costos\dfig\costosDfig.csv';

DflowD.Eolicos(2,1).nod = 6;
DflowD.Eolicos(2,1).type = 'casos\gen\dfig\Dfig_200kw_def.csv';
DflowD.Eolicos(2,1).fileG = 'casos\gen\dfig\Dfig_200kw_P_n_.csv';
DflowD.Eolicos(2,1).fileC = 'casos\costos\dfig\costosDfig.csv';

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

DflowD.Trafos(3,1).N = [-1 1];
DflowD.Trafos(3,1).TP = .1;
DflowD.Trafos(3,1).nodI = 5;
DflowD.Trafos(3,1).nodJ = 6;
DflowD.Trafos(3,1).ini = 0;
DflowD.Trafos(3,1).reg = 1;
DflowD.Trafos(3,1).cambio = .1;

% Caps
DflowD.Caps(1,1).N = [0 3];
DflowD.Caps(1,1).TP = .005;
DflowD.Caps(1,1).nod = 9;
DflowD.Caps(1,1).ini = 0;
DflowD.Caps(1,1).cambio = 1;


%% Nombres de archivos
% 
DflowD.inFilename = 'casos\PU_example\PU_example5.xls';
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
Config.outFilename = 'PU_example5_tap_cap_dfig';
Config.runNxN = false;
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

