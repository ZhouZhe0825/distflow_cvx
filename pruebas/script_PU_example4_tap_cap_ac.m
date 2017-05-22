tol = 5e-8;
CantHorasEtapa = 1;
iniEtapa = 1;

gama = .25;

gama_ini = gama;
gama_fin = gama;
TopIT = 6;
TopSecs = 60;

%% Declaracion de constantes

% Eolicos
Eol1.nod = 5;
Eol1.type = @Dfig_200kw;
Eol1.fileG = 'Dfig_200kw_P_n_.csv';
Eol1.fileC = 'costosDfig.csv';

Eolicos = [];

% Solar
Pv1.nod = 4;
Pv1.type = @PvGen_sm;
Pv1.fileG = 'pvgen.csv';
Pv1.fileC = 'costosPv.csv';

Pv2.nod = 7;
Pv2.type = @PvGen_sm;
Pv2.fileG = 'pvgen.csv';
Pv2.fileC = 'costosPv.csv';

Solares = [];

% Baterias
Bat1.nod = 5;
Bat1.type = @Bat_def;
Bat1.EIni = .5;

Baterias = [];

% Trafos
Trafo1.TP = [-2 -1 0 1 2];
Trafo1.N = .005;
Trafo1.nod = 1;
Trafo1.ini = 0;
Trafo1.cambio = 1;
Trafos = [Trafo1];

% Caps
Cap1.TP = [0 1 2 3];
Cap1.N = .005;
Cap1.nod = 9;
Cap1.ini = 0;
Cap1.cambio = 1;
Caps = [Cap1];


% Cargas
Carga1.TP = [0 1];
Carga1.pC = 0.5;
Carga1.qC = Carga1.pC *.015;
Carga1.dur = 3;
Carga1.nMultipTop = 1.1;
Carga1.nMultipLow = .75;
Carga1.nod = 5;

Carga2.TP = [0 1];
Carga2.pC = 0.25;
Carga2.qC = Carga2.pC *.015;
Carga2.dur = 2;
Carga2.nMultipTop = 1.1;
Carga2.nMultipLow = .75;
Carga2.nod = 4;

Cargas = [];

% Appliances
App1.I = 1;
App1.Pref = .8;
App1.Low = .8;
App1.Top = .8;
App1.nMultipTop = 1.05;
App1.nMultipLow = .95;
App1.alpha = 0;

App2.I = 2;
App2.Pref = .2;
App2.Low = 0;
App2.Top = .2;
App2.nMultipTop = 1.05;
App2.nMultipLow = .95;
App2.alpha = 0.1;

App = [App1, App2];

% Switches
Switches.i = [];
Switches.j = [];
Switches.cY = .1;
Switches.all = false;


tgPhi = .2;

%% Nombres de archivos
% 
outFilename_pref = 'PU_example4_tap_cap_ac';
outFilename_nxn = [outFilename_pref, '_nxn'];
outFilename_m = [outFilename_pref, '_m'];
outFilename_mat = [outFilename_pref, 'nxn2m'];

inFilename = 'PU_example4.xls';
fileCurvaCarga = 'carga_PU_example.csv';
fileUtilBetaE = '';
fileUtilBetaT = '';
fileTemp = '';
fileCostosTension = 'costosTension.csv';
fileCostosTras = 'costosTrasmision.csv';

%% Carga de datos
%% Red
[Data] = load_distflow_case(inFilename, 'bus_data_CVX', 'branch_data_CVX', Trafos, Caps, Cargas, App, Switches);

[Data] = loadCargaCuartHoraria(fileCurvaCarga, Data);


%% Generadores
% Eolicos

[Data] = cargaEolicosDefault(Data, Eolicos);

% Fotovoltaicos

[Data] = cargaPvDefault(Data, Solares);

%% Utilidad

% Configuraciones manuales

[betaE] = utilBetasE(fileUtilBetaE);

[betaT] = utilBetasT(fileUtilBetaT);

[Data] = cargaUtilDefault(Data, tgPhi, betaE, betaT, Cargas, App);

%% Parametros de Storage

% Aire Acondicionado

[Data.temp] = cargarTempInvierno(fileTemp);

[Data] = cargaACDefault(Data);

% Configuraciones manuales
Data.St.AC.I = zeros(size(Data.Red.Branch.T,1),1);
Data.St.AC.I(Data.Red.Bus.indCons) = 1;
Data.St.AC.tempLow = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * 19;
Data.St.AC.tempTop = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * 25;
Data.St.AC.tempPref = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * 22;
Data.St.AC.tempIni = Data.St.AC.I * 21;
Data.St.AC.epsilon = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * .16;
Data.St.AC.eta = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * 1666.67;
Data.St.AC.beta = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * .2;
 
% Parametros de Baterias
[Data] = cargaBatDefault(Data, Baterias);

%% Costos
[Data] = cargaCostosDefault(Data, Trafos, Caps, Switches, fileCostosTension, fileCostosTras, Solares, Eolicos);

%% Configuracion de parametros de solvers para problemas y subproblemas
% Subproblemas distribuidos
Config.iniEtapa = iniEtapa;
Config.Etapas = 4*CantHorasEtapa;
Config.workspace_var_file = outFilename_mat;

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

%% Llamada al modelo

cantTaps = length(Trafos);
cantCaps = length(Caps);
cantCargs = length(Cargas);

leyenda = ['------------------------------------ ' outFilename_mat ' ------------------------------------']

[Var_nxn, opt_nxn, DataNxN] = llamarCentralizadoNxN(Data, Config);

[Var_m, opt_m, DataM] = llamarCentralizadoM(Data, Config);

[diff_m_nxn] = checkEqualStructs(VarM2NxN(Var_m, Data), Var_nxn, 'Var_m', 'Var_nxn', 1e-5)

xlswrite([Config.workspace_var_file '_diffs.xlsx'], diff_m_nxn);

printSalidasDistflowNxN(Var_nxn, DataNxN, Config, cantTaps, cantCaps, cantCargs, outFilename_nxn, [], [], [], [], []);
printSalidasDistflowM(Var_m, DataNxN, Config, cantTaps, cantCaps, cantCargs, outFilename_m, [], [], [], [], []);

