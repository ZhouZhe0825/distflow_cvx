%% Declaracion de parametros de simulacion

% Eolicos
Eolicos = [];

% Solar
Solares = [];

% Baterias
Baterias = [];

% Trafos
Trafos = [];

% Caps
Caps = [];


% Cargas
Cargas = [];

% Appliances
App(1).I = 1;
App(1).Pref = .8;
App(1).Low = .8;
App(1).Top = .8;
App(1).nMultipTop = 1.05;
App(1).nMultipLow = .95;
App(1).alpha = 0;
App(1).tgPhi = .2;

App(2).I = 2;
App(2).Pref = .2;
App(2).Low = 0;
App(2).Top = .2;
App(2).nMultipTop = 1.05;
App(2).nMultipLow = .95;
App(2).alpha = 0.1;
App(2).tgPhi = .2;

% Aires acondicionados
ACs = [];

% Switches
Switches.i = [];
Switches.j = [];
Switches.cY = .1;
Switches.all = false;


%% Nombres de archivos
% 
inFilename = 'casos\PU_example\PU_example4.xls';
fileCurvaCarga = 'casos\PU_example\cargas\carga_PU_example.csv';
fileUtilBetaT = 'casos\util\betaT.csv';
fileTemp = 'casos\temp\tempInvierno.csv';
fileCostosTension = 'casos\costos\tension\costosTension.csv';
fileCostosTras = 'casos\costos\trasmision\costosTrasmision.csv';

%% Configuracion de simulacion
Config.iniEtapa = 1;
Config.Etapas = 1;
Config.outFilename = 'PU_example4';
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
