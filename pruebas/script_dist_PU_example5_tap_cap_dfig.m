%% Declaracion de parametros de simulacion

DflowD = cargaDflowDDefault();

% Eolicos
Eol1.nod = 5;
Eol1.type = 'casos\gen\dfig\Dfig_200kw_def.csv';
Eol1.fileG = 'casos\gen\dfig\Dfig_200kw_P_n_.csv';
Eol1.fileC = 'casos\costos\dfig\costosDfig.csv';

Eol2.nod = 6;
Eol2.type = 'casos\gen\dfig\Dfig_200kw_def.csv';
Eol2.fileG = 'casos\gen\dfig\Dfig_200kw_P_n_.csv';
Eol2.fileC = 'casos\costos\dfig\costosDfig.csv';

DflowD.Eolicos = [Eol1;Eol2];

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
DflowD.inFilename = 'casos\PU_example\PU_example5.xls';
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
Config.outFilename = 'PU_example5_tap_cap_dfig_dist';
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

