cambioTap = 1;
cambioCap = 1;
minpC = 0;
minqC = minpC * .15;
minr = 0;
tol = 5e-8;
delta = 0.004;
m = .01;
CantHorasEtapa = .25;
iniEtapa = 1;

gama = .25;

gama_ini = gama;
gama_fin = gama;
TopIT = 6;
TopSecs = 60;

%% Declaracion de constantes

NodosGeneracionEolica = [];
NodosGeneracionSolar = [];
NodosBaterias = [];

% Trafos
Trafo1.TP = [0];
Trafo1.N = .005;
Trafo1.nod = 1;
Trafo1.ini = 0;
Trafos = [Trafo1];

% Caps
Cap1.TP = [0];
Cap1.N = .005;
Cap1.nod = 9;
Cap1.ini = 0;
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
Switches.all = false;


tgPhi = .2;

rhoQ_ct = 1.5;
rhoP_ct = 1;

cv_ct = .1;
cr_ct = .1;
cY = .1;
EIni_ct = .5;
betaT = .2;

%% Nombres de archivos
% 
outFilename_pref = 'PU_example5';
outFilename_nxn = [outFilename_pref, '_nxn'];
outFilename_m = [outFilename_pref, '_m'];
outFilename_mat = [outFilename_pref, 'nxn2m'];

inFilename = 'PU_example5.xls';
fileCurvaCarga = 'carga_PU_example.xlsx';


%% Carga de datos
%% Red
[Data] = load_distflow_case(inFilename, 'bus_data_CVX', 'branch_data_CVX', Trafos, Caps, Cargas, App, Switches, cambioTap, cambioCap, cY);

[Data] = loadCargaCuartHoraria(fileCurvaCarga, Data, 'Ppu', 'Qpu', minpC, minqC);


%% Generadores
% Eolicos
[Data] = Dfig_200kw(Data, NodosGeneracionEolica);

[vVel, Data.temp] = cargarVelocidadVientoInvierno();

[Data.Gen.DFIG.P_mec, Data.Gen.DFIG.n_]= calculoPotenciaEolica(vVel, ...
	Data.Gen.DFIG.vmm, Data.Gen.DFIG.vm, Data.Gen.DFIG.vM, Data.Gen.DFIG.vMM, ...
	Data.Gen.DFIG.Omega, Data.Gen.DFIG.G, Data.Gen.DFIG.P_nMec, Data.Gen.DFIG.Np, ...
	Data.Gen.DFIG.R_, Data.Gen.DFIG.rho, Data.Gen.DFIG.ws, Data.Gen.DFIG.c1, ...
	Data.Gen.DFIG.c2, Data.Gen.DFIG.c3, Data.Gen.DFIG.c4, Data.Gen.DFIG.c5, ...
	Data.Gen.DFIG.c6, Data.Gen.DFIG.c7, Data.Gen.DFIG.lambda_opt);

% Fotovoltaicos

[Data] = PvGen_sm(Data,NodosGeneracionSolar);

%% Utilidad

% Configuraciones manuales

utilCarg = ones(1,96);
utilCarg(1) = .125;
utilCarg(2) = .375;
utilCarg(3) = .625;
utilCarg(4) = .875;

[Data] = cargarUtilDefault(Data, tgPhi, betaT, utilCarg, Cargas, App);

%% Parametros de Storage

% Configuraciones manuales

% Aire Acondicionado
Data.dt = .25;
Data.St.AC.I = zeros(size(Data.Red.Branch.T,1),1);
Data.St.AC.tempLow = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * 19;
Data.St.AC.tempTop = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * 25;
Data.St.AC.tempPref = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * 22;
Data.St.AC.tempIni = Data.St.AC.I * 21;
Data.St.AC.epsilon = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * .16;
Data.St.AC.eta = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * 1666.67;
Data.St.AC.beta = repmat(Data.St.AC.I, [1, size(Data.temp,2)]) * .2;
 
% Parametros de Baterias
[Data] = Bat_def(Data,EIni_ct,NodosBaterias);

%% Costos
[mHor, piPTrasHor] = mRhoHorario();

[Data] = cargaCostosDefault(Data, mHor, piPTrasHor, rhoP_ct, rhoQ_ct, m, delta);

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

[diff_m_nxn] = checkEqualStructs(Var_m, Var_nxn, 'Var_m', 'Var_nxn', 1e-5)

xlswrite([Config.workspace_var_file '_diffs.xlsx'], diff_m_nxn);

% printSalidasDistflow(Var_nxn, DataNxN, Config, cantTaps, cantCaps, cantCargs, outFilename_nxn, [], [], [], [], []);
% printSalidasDistflow(Var_m, DataNxN, Config, cantTaps, cantCaps, cantCargs, outFilename_m, [], [], [], [], []);

