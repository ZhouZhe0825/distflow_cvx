%% Declaracion de parametros de simulacion

DflowD = cargaDflowDDefault();

% Solar
DflowD.Solares(1,1).nod = 4;
DflowD.Solares(1,1).type = 'casos\gen\pv\pv_def.csv';
DflowD.Solares(1,1).fileG = 'casos\gen\pv\pvgen.csv';
DflowD.Solares(1,1).fileC = 'casos\costos\pv\costosPv.csv';

% Baterias
DflowD.Baterias(1,1).nod = 5;
DflowD.Baterias(1,1).type = 'casos\bat\bat_def.csv';
DflowD.Baterias(1,1).EIni = .5;

% Eolicos
DflowD.Eolicos(1,1).nod = 5;
DflowD.Eolicos(1,1).type = 'casos\gen\dfig\Dfig_1mw_def.csv';
DflowD.Eolicos(1,1).fileG = 'casos\gen\dfig\Dfig_1mw_P_n_.csv';
DflowD.Eolicos(1,1).fileC = 'casos\costos\dfig\costosDfig_1.csv';

% Generador Basico
DflowD.GenBas(1,1).nod = 7;
DflowD.GenBas(1,1).fileG = 'casos\gen\basic\genbas.csv';
DflowD.GenBas(1,1).fileC = 'casos\costos\basic\costosBasic.csv';

% Switches
DflowD.Switches.i = [5 7];
DflowD.Switches.j = [6 8];
DflowD.Switches.cY = 1;
DflowD.Switches.all = false;

% Trafos
DflowD.Trafos(1,1).N = [2 8];
DflowD.Trafos(1,1).TP = .025;
DflowD.Trafos(1,1).nodI = 1;
DflowD.Trafos(1,1).nodJ = 2;
DflowD.Trafos(1,1).ini = 2;
DflowD.Trafos(1,1).reg = 0;
DflowD.Trafos(1,1).cambio = 1;

DflowD.Trafos(2,1).N = [-1 1];
DflowD.Trafos(2,1).TP = .1;
DflowD.Trafos(2,1).nodI = 4;
DflowD.Trafos(2,1).nodJ = 5;
DflowD.Trafos(2,1).ini = 1;
DflowD.Trafos(2,1).reg = 1;
DflowD.Trafos(2,1).cambio = .1;

% Caps
DflowD.Caps(1,1).N = [0 3];
DflowD.Caps(1,1).TP = .005;
DflowD.Caps(1,1).nod = 9;
DflowD.Caps(1,1).ini = 3;
DflowD.Caps(1,1).cambio = 1;

% Aires acondicionados
DflowD.ACs(1,1).fileT = 'casos\PU_example\ac\ac.csv';
DflowD.ACs(1,1).tempIni = 21;
DflowD.ACs(1,1).epsilon = .16;
DflowD.ACs(1,1).eta = 1666.67;

% Cargas
DflowD.Cargas(1,1).pC = 0.5;
DflowD.Cargas(1,1).qC = DflowD.Cargas(1,1).pC *.015;
DflowD.Cargas(1,1).dur = 3;
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
DflowD.Cargas(2,1).nod = 4;

%% Nombres de archivos
% 
DflowD.inFilename = 'casos\PU_example\PU_example4.xls';
DflowD.fileCurvaCarga = 'casos\PU_example\cargas\carga_PU_example.csv';
DflowD.fileUtilBetaT = 'casos\util\betaT.csv';
DflowD.utilOptFuncCuad = false;
DflowD.fileTemp = 'casos\temp\tempInvierno.csv';
DflowD.fileCostosTension = 'casos\costos\tension\costosTension.csv';
DflowD.fileCostosTras = 'casos\costos\trasmision\costosTrasmision.csv';

%% Configuracion de simulacion
iniEtapa = 1;
CantHorasEtapa = 1;

Config = [];
Config.iniEtapa = iniEtapa;
Config.Etapas = 4*CantHorasEtapa;
Config.outFilename = 'PU_example4_all_dist';
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

