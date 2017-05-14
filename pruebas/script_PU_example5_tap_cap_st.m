tol = 5e-8;
CantHorasEtapa = 1;
iniEtapa = 1;

gama = .25;

gama_ini = gama;
gama_fin = gama;
TopIT = 6;
TopSecs = 60;

%% Declaracion de constantes

NodosGeneracionEolica = [];
NodosGeneracionSolar = [];
NodosBaterias = [5];

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
Carga1.dur = 4;
Carga1.nMultipTop = 1.1;
Carga1.nMultipLow = .75;
Carga1.nod = 5;

Carga2.TP = [0 1];
Carga2.pC = 0.25;
Carga2.qC = Carga2.pC *.015;
Carga2.dur = 2;
Carga2.nMultipTop = 1.1;
Carga2.nMultipLow = .75;
Carga2.nod = 9;

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
EIni_ct = .5;

%% Nombres de archivos
% 
outFilename_pref = 'PU_example5_tap_cap_st';
outFilename_nxn = [outFilename_pref, '_nxn'];
outFilename_m = [outFilename_pref, '_m'];
outFilename_mat = [outFilename_pref, 'nxn2m'];

inFilename = 'PU_example5.xls';
fileCurvaCarga = 'carga_PU_example.xlsx';
fileP_mec = '';
filePPvg = '';
fileUtilBeta = '';
fileTemp = '';
fileCostosTension = '';
fileCostosPv = '';
fileCostosDfig = '';
fileCostosTras = '';

%% Carga de datos
%% Red
[Data] = load_distflow_case(inFilename, 'bus_data_CVX', 'branch_data_CVX', Trafos, Caps, Cargas, App, Switches);

[Data] = loadCargaCuartHoraria(fileCurvaCarga, Data, 'Ppu', 'Qpu');


%% Generadores
% Eolicos

[Data] = cargaEolicosDefault(Data);

[P_mec, n_] = p_mecN_(fileP_mec);

[Data] = Dfig_200kw(Data, NodosGeneracionEolica, P_mec, n_);

% Fotovoltaicos

[Data] = cargaPvDefault(Data);

[pPvg] = pPvgs(filePPvg);

[Data] = PvGen_sm(Data,NodosGeneracionSolar, pPvg);

%% Utilidad

% Configuraciones manuales

[utilCarg, betaT] = utilBetas(fileUtilBeta);

[Data] = cargaUtilDefault(Data, tgPhi, utilCarg, betaT, Cargas, App);

%% Parametros de Storage

% Aire Acondicionado

[Data.temp] = cargarTempInvierno(fileTemp);

[Data] = cargaACDefault(Data);

% Configuraciones manuales
 
% Parametros de Baterias
[Data] = Bat_def(Data,EIni_ct,NodosBaterias);

%% Costos
[mHor, delta, cdv] = costosTension(fileCostosTension);

[piPTrasHor, piQmtrasHor, piQMtrasHor] = costosTrasmision(fileCostosTras);

[rhopPv,rhomqPv,rhoMqPv] = costosPv(fileCostosPv);

[rhopWi,rhomqWi,rhoMqWi] = costosDfig(fileCostosDfig);

[Data] = cargaCostosDefault(Data, Trafos, Caps, Switches, mHor, cdv, delta, rhopPv, rhomqPv, rhoMqPv, rhopWi, rhomqWi, rhoMqWi, piPTrasHor, piQmtrasHor, piQMtrasHor);

%% Configuracion de parametros de solvers para problemas y subproblemas
% Subproblemas distribuidos
Config.iniEtapa = iniEtapa;
Config.Etapas = 4*CantHorasEtapa;
Config.workspace_var_file = outFilename_mat;

% Centralizado
Config.Centr{1, 1} = 'MSK_DPAR_MIO_TOL_REL_RELAX_INT';
Config.Centr{1, 2} = 0.01;
Config.Centr{2, 1} = 'MSK_DPAR_OPTIMIZER_MAX_TIME';
Config.Centr{2, 2} = 600;

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

