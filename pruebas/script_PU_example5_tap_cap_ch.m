%% Declaracion de parametros de simulacion

% Eolicos
Eol1.nod = 5;
Eol1.type = 'casos\gen\dfig\Dfig_200kw_def.csv';
Eol1.fileG = 'casos\gen\dfig\Dfig_200kw_P_n_.csv';
Eol1.fileC = 'casos\costos\dfig\costosDfig.csv';

Eol2.nod = 6;
Eol2.type = 'casos\gen\dfig\Dfig_200kw_def.csv';
Eol2.fileG = 'casos\gen\dfig\Dfig_200kw_P_n_.csv';
Eol2.fileC = 'casos\costos\dfig\costosDfig.csv';

Eolicos = [];

% Solar
Pv1.nod = 4;
Pv1.type = 'casos\gen\pv\pv_def.csv';
Pv1.fileG = 'casos\gen\pv\pvgen.csv';
Pv1.fileC = 'casos\costos\pv\costosPv.csv';

Pv2.nod = 7;
Pv2.type = 'casos\gen\pv\pv_def.csv';
Pv2.fileG = 'casos\gen\pv\pvgen.csv';
Pv2.fileC = 'casos\costos\pv\costosPv.csv';

Solares = [];

% Baterias
Bat1.nod = 5;
Bat1.type = 'casos\bat\bat_def.csv';
Bat1.EIni = .5;

Baterias = [];

% Trafos
Trafo1.N = [-2 2];
Trafo1.TP = .005;
Trafo1.nod = 1;
Trafo1.ini = 0;
Trafo1.cambio = 1;
Trafos = [Trafo1];

% Caps
Cap1.N = [0 3];
Cap1.TP = .005;
Cap1.nod = 9;
Cap1.ini = 0;
Cap1.cambio = 1;
Caps = [Cap1];


% Cargas
Carga1.TP = [0 1];
Carga1.pC = 0.5;
Carga1.qC = Carga1.pC *.015;
Carga1.dur = 4;
Carga1.nMultipTop = 1.1;
Carga1.nMultipLow = .75;
Carga1.fileU = 'casos\util\betaE.csv';
Carga1.nod = 5;

Carga2.TP = [0 1];
Carga2.pC = 0.25;
Carga2.qC = Carga2.pC *.015;
Carga2.dur = 2;
Carga2.nMultipTop = 1.1;
Carga2.nMultipLow = .75;
Carga2.fileU = 'casos\util\betaE.csv';
Carga2.nod = 9;

Cargas = [Carga2];

% Appliances
App1.I = 1;
App1.Pref = .8;
App1.Low = .8;
App1.Top = .8;
App1.nMultipTop = 1.05;
App1.nMultipLow = .95;
App1.alpha = 0;
App1.tgPhi = .2;

App2.I = 2;
App2.Pref = .2;
App2.Low = 0;
App2.Top = .2;
App2.nMultipTop = 1.05;
App2.nMultipLow = .95;
App2.alpha = 0.1;
App2.tgPhi = .2;

App = [App1;App2];

% Aires acondicionados
Ac1.nod = [2;3;4;5;6;7;8;9];
Ac1.fileT = 'casos\PU_example\ac\ac.csv';
Ac1.tempIni = 21;
Ac1.epsilon = .16;
Ac1.eta = 1666.67;

ACs = [];

% Switches
Switches.i = [];
Switches.j = [];
Switches.cY = .1;
Switches.all = false;


%% Nombres de archivos
% 
inFilename = 'casos\PU_example\PU_example5.xls';
fileCurvaCarga = 'casos\PU_example\cargas\carga_PU_example.csv';
fileUtilBetaT = 'casos\util\betaT.csv';
fileTemp = 'casos\temp\tempInvierno.csv';
fileCostosTension = 'casos\costos\tension\costosTension.csv';
fileCostosTras = 'casos\costos\trasmision\costosTrasmision.csv';

%% Configuracion de simulacion
iniEtapa = 1;
CantHorasEtapa = 1;

Config.iniEtapa = iniEtapa;
Config.Etapas = 4*CantHorasEtapa;
Config.outFilename = 'PU_example5_tap_cap_ch';
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

[Data] = loadData(Trafos, Caps, Switches, App, Eolicos, Solares, Baterias, Cargas, ACs, inFilename, fileCurvaCarga, fileUtilBetaT, fileTemp, fileCostosTension, fileCostosTras);

%% Correr Simulacion

runSimulation(Data, Trafos, Caps, Switches, App, Eolicos, Solares, Baterias, Cargas, ACs, inFilename, fileCurvaCarga, fileUtilBetaT, fileTemp, fileCostosTension, fileCostosTras, Config)

