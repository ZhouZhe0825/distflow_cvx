alphaMult = 0;
cambioRed = 1.5;
cambioTap = 1;
cambioCap = 1;
cg1pc = 0.5;
cg1durHor = 4;
cg2pc = 0.25;
cg2durHor = 2;
util = true;
minpC = 0;
minr = 0;
tol = 5e-8;
delta = 0.004;
m = .01;
CantHorasEtapa = .5;
iniEtapa = 1; %7:00

factorReducEol = .25;

gama_ini = .25;
gama_fin = .25;
TopIT = 6;
TopSecs = 60;
util=true;

%% Declaracion de constantes

% NodosGeneracionEolica = [5];
NodosGeneracionEolica = [];

% Trafos
Trafo1.TP = [3];
Trafo1.N = .00005;
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
Carga1.pC = cg1pc;
Carga1.qC = cg1pc*.015;
Carga1.dur = cg1durHor;
Carga1.nod = 5;

Carga2.TP = [0 1];
Carga2.pC = cg2pc;
Carga2.qC = cg2pc*.015;
Carga2.dur = cg2durHor;
Carga2.nod = 9;

% Cargas = [Carga1;Carga2];
% Cargas = [Carga2];
Cargas = [];


multipCarga = ones(9, 96);
rhoQ_ct = 0.1;
rhoP_ct = 0.1;

dual = 1;



cv_ct = .1;
cr_ct = .1;
lTop_ct = 5;

uLow_ct = .95;
uTop_ct = 1.05;

iniEstado = 1;
%% Nombres de archivos
% 
outFilename_pref = 'PU_example4_St';
outFilename_c = [outFilename_pref, '_nxn'];
outFilename_r = [outFilename_pref, '_m'];
outFilename_mat = [outFilename_pref, 'nxn2m.mat'];

% outFilenamePre = ['.\paper\prueba_inv_'];
% outFilenameMid = ['d_g_25_h0_1_DFig_Sto_Ch_betLess'];
% outFilenameSuf = [ ...
%     '_tol' num2str(tol) '_it' num2str(TopIT)];

% outFilename_w_sw_dist = [outFilenamePre '_d_' outFilenameMid outFilenameSuf ''];

inFilename = 'PU_example4.xls';
% fileCurvaCarga = 'carga_subredLP_trafo_red_93-73_no_rnd.xlsx';

%% Carga de datos
%% Red
[Data] = load_distflow_case(inFilename, 'bus_data_CVX', 'branch_data_CVX');

% Se asume que hay hasta un Trafo por nodo
Data.Red.Bus.TapLow = Data.Red.Bus.uLow * 0;
Data.Red.Bus.TapTop = Data.Red.Bus.TapLow;
Data.Red.Bus.TapIni = Data.Red.Bus.TapLow;
Data.Red.Bus.indTap = Data.Red.Bus.TapLow;

for estT = 1:size(Trafos,1)
	Data.Red.Bus.Ntr(Trafos(estT).nod) = Trafos(estT).N;
	Data.Red.Bus.TapLow(Trafos(estT).nod) = min(Trafos(estT).TP);
	Data.Red.Bus.TapTop(Trafos(estT).nod) = max(Trafos(estT).TP);
	Data.Red.Bus.TapIni(Trafos(estT).nod) = Trafos(estT).ini;
    Data.Red.Bus.indTap(Trafos(estT).nod) = 1;
end

% Se asume que hay hasta un Capacitor por nodo
Data.Red.Bus.CapLow = Data.Red.Bus.uLow * 0;
Data.Red.Bus.CapTop = Data.Red.Bus.CapLow;
Data.Red.Bus.CapIni = Data.Red.Bus.CapLow;
Data.Red.Bus.indCap = Data.Red.Bus.CapLow;

for estC = 1:size(Caps,1)
	Data.Red.Bus.CapLow(Caps(estC).nod) = min(Caps(estC).TP);
	Data.Red.Bus.CapTop(Caps(estC).nod) = max(Caps(estC).TP);
	Data.Red.Bus.Ncp(Caps(estC).nod) = Caps(estC).N;
	Data.Red.Bus.CapIni(Caps(estC).nod) = Caps(estC).ini;
    Data.Red.Bus.indCap(Caps(estC).nod) = 1;
end

Data.Red.cambioRed = cambioRed;
Data.Red.cambioTap = cambioTap;
Data.Red.cambioCap = cambioCap;

% Configuraciones manuales
Data.Red.Bus.pCLow = repmat(Data.Red.Bus.pCLow, 1, 96);
Data.Red.Bus.qCLow = repmat(Data.Red.Bus.qCLow, 1, 96);
Data.Red.Bus.pCLow = Data.Red.Bus.pCLow.*multipCarga;
Data.Red.Bus.qCLow = Data.Red.Bus.qCLow.*multipCarga;

Data.Red.Bus.alpha = repmat(full(Data.Red.Bus.alpha), [1,2]);
Data.Red.Bus.alpha(:,1) = 0;
Data.Red.Bus.alpha(:,2) = .1;

Data.Red.Bus.indCons = intersect( find(Data.Red.Bus.pCLow(:,1) ~= 0), find(Data.Red.Bus.qCLow(:,1) ~= 0));

auxV = Data.Red.Bus.uTop(Data.Red.Bus.v0);
Data.Red.Bus.uLow = ones(size(Data.Red.Bus.uLow))*uLow_ct;
Data.Red.Bus.uTop = ones(size(Data.Red.Bus.uLow))*uTop_ct;
Data.Red.Bus.uLow(Data.Red.Bus.v0) = auxV;
Data.Red.Bus.uTop(Data.Red.Bus.v0) = auxV;

Data.Red.Branch.lTop(:,:) = lTop_ct;

Data.Red.Branch.yTop = Data.Red.Branch.T;
Data.Red.Branch.yLow = Data.Red.Branch.T;

Data.Red.Branch.yLow(8, 9) = 0;
Data.Red.Branch.yLow(9, 8) = 0;
    
%% Clientes no interrumpibles
Data.ClNI.pC = Data.Red.Bus.uLow * 0;
Data.ClNI.qC = Data.ClNI.pC;
Data.ClNI.d = Data.ClNI.pC;
Data.ClNI.I = Data.ClNI.pC;
for estCg = 1:size(Cargas,1)
	Data.ClNI.pC(Cargas(estCg).nod) = Cargas(estCg).pC;
	Data.ClNI.qC(Cargas(estCg).nod) = Cargas(estCg).qC;
	Data.ClNI.d(Cargas(estCg).nod) = Cargas(estCg).dur;
    Data.ClNI.I(Cargas(estCg).nod) = 1;
end
	
%% Generadores
% Eolicos
[Data] = cargaDatosEolicos(Data, NodosGeneracionEolica, 0);

[vVel, temp] = cargarVelocidadVientoInvierno();

[Data.Gen.DFIG.P_mec, Data.Gen.DFIG.n_]= calculoPotenciaEolica(vVel, ...
	Data.Gen.DFIG.vmm, Data.Gen.DFIG.vm, Data.Gen.DFIG.vM, Data.Gen.DFIG.vMM, ...
	Data.Gen.DFIG.Omega, Data.Gen.DFIG.G, Data.Gen.DFIG.P_nMec, Data.Gen.DFIG.Np, ...
	Data.Gen.DFIG.R_, Data.Gen.DFIG.rho, Data.Gen.DFIG.ws, Data.Gen.DFIG.c1, ...
	Data.Gen.DFIG.c2, Data.Gen.DFIG.c3, Data.Gen.DFIG.c4, Data.Gen.DFIG.c5, ...
	Data.Gen.DFIG.c6, Data.Gen.DFIG.c7, Data.Gen.DFIG.lambda_opt);
Data.temp = temp;

% Fotovoltaicos

% Trasmision
% Data.Red.Bus.Q0Top = .0955;
% Data.Red.Bus.Q0Low = 0;
% Data.Red.Bus.P0Top = .416;
% Data.Red.Bus.P0Low = 0;

Data.Gen.Tras.pgLow = Data.Red.Bus.pCLow * 0;
Data.Gen.Tras.qgLow = Data.Gen.Tras.pgLow;
Data.Gen.Tras.pgTop = Data.Gen.Tras.pgLow;
Data.Gen.Tras.qgTop = Data.Gen.Tras.pgLow;

Data.Gen.Tras.pgLow(Data.Red.Bus.v0,:) = Data.Red.Bus.P0Low;
Data.Gen.Tras.qgLow(Data.Red.Bus.v0,:) = Data.Red.Bus.Q0Low;
Data.Gen.Tras.pgTop(Data.Red.Bus.v0,:) = Data.Red.Bus.P0Top;
Data.Gen.Tras.qgTop(Data.Red.Bus.v0,:) = Data.Red.Bus.Q0Top;

Data.Gen.Tras.I = zeros(size(Data.Red.Branch.T,1),1);
Data.Gen.Tras.I(Data.Red.Bus.v0) = 1;


% Configuraciones manuales
Data.Gen.DFIG.P_nMec = 150000;
Data.Gen.DFIG.cv = Data.Gen.DFIG.cv*1;
Data.Gen.DFIG.cr = Data.Gen.DFIG.cr*1;
temp = ones(size(temp))*temp(1);

vVel = ones(size(vVel))*15;

[Data.Gen.DFIG.P_mec, Data.Gen.DFIG.n_]= calculoPotenciaEolica(vVel, ...
	Data.Gen.DFIG.vmm, Data.Gen.DFIG.vm, Data.Gen.DFIG.vM, Data.Gen.DFIG.vMM, ...
	Data.Gen.DFIG.Omega, Data.Gen.DFIG.G, Data.Gen.DFIG.P_nMec, Data.Gen.DFIG.Np, ...
	Data.Gen.DFIG.R_, Data.Gen.DFIG.rho, Data.Gen.DFIG.ws, Data.Gen.DFIG.c1, ...
	Data.Gen.DFIG.c2, Data.Gen.DFIG.c3, Data.Gen.DFIG.c4, Data.Gen.DFIG.c5, ...
	Data.Gen.DFIG.c6, Data.Gen.DFIG.c7, Data.Gen.DFIG.lambda_opt);

Data.Gen.DFIG.P_mec = Data.Gen.DFIG.P_mec';
Data.Gen.DFIG.n_ = Data.Gen.DFIG.n_';
Data.Gen.DFIG.P_mec = Data.Gen.DFIG.P_mec/2;

% Data.Gen.Pv.I(7) = 1;
Data.Gen.Pv.pPvg = ones(length(Data.Red.Branch.T),1).*Data.Gen.Pv.I*100; %potencia de generacion del solar, por ahora constante

indSn = find(Data.Gen.Pv.I == 1);
Data.Gen.Pv.sTop(indSn) = 0.8;											
Data.Gen.Pv.pgTop(indSn) = 0.0075;											
Data.Gen.Pv.cv(indSn) = cv_ct;
Data.Gen.Pv.cr(indSn) = cr_ct;
Data.temp = temp;



%% Utilidad

% Configuraciones manuales

utilCarg = utilidadCarga();
if util
    Data.Util.Func = 2;

    Data.Util.betaT = zeros(size(Data.Red.Bus.pCLow,1),1);
    Data.Util.betaT(Data.Red.Bus.indCons) = .2;
    
    Data.Util.betaE = zeros(size(Data.Red.Bus.pCLow));
    for estCg = 1:size(Cargas,1)
        Data.Util.betaE(Cargas(estCg).nod,:) = utilCarg';
    end
    
    Data.Util.aT = ones(size(Data.Red.Bus.pCLow,1),1)*0;
    Data.Util.aE = ones(size(Data.Red.Bus.pCLow,1),1)*0;
    
    
    Data.Util.nMultipTop = 1.05;
    Data.Util.nMultipLow = .95;
    Data.Util.tgPhi = .2;
    
else
    Data.Util.Func = 0;
end

if Data.Util.Func ~= 0
	Data.Util.pzCnPref = repmat(full(Data.Red.Bus.pCLow), [1,1,2]);		
	Data.Util.pzCnLow = repmat(full(Data.Red.Bus.pCLow), [1,1,2]) * Data.Util.nMultipLow;
	Data.Util.pzCnTop = repmat(full(Data.Red.Bus.pCLow), [1,1,2]) * Data.Util.nMultipTop;

	Data.Util.pzCnPref(:,:,1) = Data.Util.pzCnPref(:,:,1) * .8;		
	Data.Util.pzCnLow(:,:,1) = Data.Util.pzCnLow(:,:,1) * .8;
	Data.Util.pzCnTop(:,:,1) = Data.Util.pzCnTop(:,:,1) * .8;

	Data.Util.pzCnPref(:,:,2) = Data.Util.pzCnPref(:,:,2) * .2;		
	Data.Util.pzCnLow(:,:,2) = Data.Util.pzCnLow(:,:,2) * 0;
	Data.Util.pzCnTop(:,:,2) = Data.Util.pzCnTop(:,:,2) * .2;

	Data.Util.pzCnPrefE = full(Data.Red.Bus.pCLow) * 0;
	Data.Util.pzCnLowE = full(Data.Red.Bus.pCLow) * 0;
	Data.Util.pzCnTopE = full(Data.Red.Bus.pCLow) * 0;

	Data.Util.qzCnPrefE = full(Data.Red.Bus.pCLow) * 0;
	Data.Util.qzCnLowE = full(Data.Red.Bus.pCLow) * 0;
	Data.Util.qzCnTopE = full(Data.Red.Bus.pCLow) * 0;

	Data.Util.pzCnPrefE = Data.ClNI.pC;
	Data.Util.pzCnLowE = Data.ClNI.pC*.75;
	Data.Util.pzCnTopE = Data.ClNI.pC*1.1;

	Data.Util.qzCnPrefE = Data.ClNI.qC;
	Data.Util.qzCnLowE = Data.ClNI.qC*.75;
	Data.Util.qzCnTopE = Data.ClNI.qC*1.1;
end

%% Parametros de Storage

% Configuraciones manuales

% Aire Acondicionado
Data.St.AC.tempLow = zeros(size(Data.Red.Branch.T,1),1); % temperatura minima, por nodo
Data.St.AC.tempLow(Data.Red.Bus.indCons) = 19;
Data.St.AC.tempTop = zeros(size(Data.Red.Branch.T,1),1); % temperatura maxima, por nodo
Data.St.AC.tempTop(Data.Red.Bus.indCons) = 25;
Data.St.AC.tempPref = zeros(size(Data.Red.Branch.T,1),1); % temperatura preferida, por nodo
Data.St.AC.tempPref(Data.Red.Bus.indCons) = 22;
Data.St.AC.tempIni =  zeros(size(Data.Red.Branch.T,1),1); % temperatura preferida, por nodo
Data.St.AC.tempIni(Data.Red.Bus.indCons) = 21;
Data.St.AC.epsilon = zeros(size(Data.Red.Branch.T,1),1); % epsilon por nodo
Data.St.AC.epsilon(Data.Red.Bus.indCons) = .16;
% dt = 15 * 60;
Data.dt = .25;
Data.St.AC.eta = zeros(size(Data.Red.Branch.T,1),1); % eta por nodo
Data.St.AC.eta(Data.Red.Bus.indCons) = 1666.67;
 
Data.St.AC.beta = zeros(size(Data.Red.Branch.T,1),1); %betaAC por nodo
Data.St.AC.beta(Data.Red.Bus.indCons) = 0;
Data.St.AC.a = zeros(size(Data.Red.Branch.T,1),1); %aAC por nodo
 
% Parametros de Baterias
Data.St.Bat.I = zeros(size(Data.Red.Branch.T,1),1);
% Data.St.Bat.I(5) = 1;
Data.St.Bat.cv = Data.St.Bat.I * .001; % temperatura minima, por nodo
Data.St.Bat.cr = Data.St.Bat.I * .001; % temperatura maxima, por nodo

Data.St.Bat.epsilon = Data.St.Bat.I * .1; % epsilon por nodo
Data.St.Bat.eta = Data.St.Bat.I; % eta por nodo

Data.St.Bat.pgTop = Data.St.Bat.I * .5; % temperatura minima, por nodo
Data.St.Bat.pgLow = Data.St.Bat.I * -.5; % temperatura maxima, por nodo
Data.St.Bat.sTop = Data.St.Bat.I * .1; % temperatura minima, por nodo

Data.St.Bat.ETop = Data.St.Bat.I * 2; % temperatura minima, por nodo
Data.St.Bat.EIni = Data.St.Bat.I * .5; % temperatura minima, por nodo
Data.St.Bat.ELow = Data.St.Bat.I * 0; % temperatura maxima, por nodo

Data.St.Bat.beta = Data.St.Bat.I * .2; %betaAC por nodo
Data.St.Bat.wU = Data.St.Bat.I; %betaAC por nodo
Data.St.Bat.wOm = Data.St.Bat.I; %aAC por nodo

Data.St.Bat.m1 = Data.St.Bat.I * 1; %aAC por nodo
Data.St.Bat.m2 = Data.St.Bat.I * .75; %aAC por nodo
Data.St.Bat.m3 = Data.St.Bat.I * .5; %aAC por nodo

Data.St.Bat.kapa = .2; %aAC por nodo
Data.St.Bat.gama = .6; %aAC por nodo


%% Costos
[Data] = cargaCostosDefault(Data);

% Configuraciones manuales
[mHor, piPTrasHor] = mRhoHorario();

Data.Cost.piPTras = Data.Cost.piPTras * rhoP_ct * piPTrasHor';
Data.Cost.piQmtras = Data.Cost.piQmtras * rhoQ_ct * piPTrasHor';
Data.Cost.piQMtras = Data.Cost.piQMtras * rhoQ_ct * piPTrasHor';

Data.Cost.m = m*mHor';
Data.Cost.delta = delta;

%% Configuracion de parametros de solvers para problemas y subproblemas
% problema Centralizado e inicial hallados

Var_centr = [];
opt_centr = [];
Var_ini = [];
opt_ini = [];

% if exist([outFilename_w_sw_dist, '_var.mat'], 'file')
%     load([outFilename_w_sw_dist, '_var.mat'], 'Var_centr', 'opt_centr', 'Var_ini', 'opt_ini');
% end


% Subproblemas distribuidos
Config.iniEtapa = iniEtapa;
Config.Etapas = 4*CantHorasEtapa;
Config.workspace_var_file = outFilename_mat;

Config.Distr.gama_ini = gama_ini;
Config.Distr.gama_fin = gama_fin;
Config.Distr.TopIT = TopIT;
Config.Distr.tol = tol;

Config.SubP{1,1} = 'MSK_DPAR_MIO_MAX_TIME';
Config.SubP{1,2} = TopSecs/2;
Config.SubP{2,1} = 'MSK_DPAR_OPTIMIZER_MAX_TIME';
Config.SubP{2,2} = TopSecs;
Config.SubP{3,1} = 'MSK_DPAR_MIO_TOL_REL_RELAX_INT';
Config.SubP{3,2} = 0.01;
Config.SubP{4,1} = 'MSK_DPAR_MIO_DISABLE_TERM_TIME';
Config.SubP{4,2} = 10;
Config.SubP{5,1} = 'MSK_DPAR_MIO_NEAR_TOL_REL_GAP';
Config.SubP{5,2} = 0.025;
Config.SubP{6,1} = 'MSK_IPAR_OPF_WRITE_PROBLEM';
Config.SubP{6,2} = 'MSK_ON';
Config.SubP{7,1} = 'MSK_IPAR_OPF_WRITE_HEADER';
Config.SubP{7,2} = 'MSK_ON';
Config.SubP{8,1} = 'MSK_IPAR_INFEAS_REPORT_AUTO';
Config.SubP{8,2} = 'MSK_ON';
Config.SubP{9,1} = 'MSK_SPAR_PARAM_WRITE_FILE_NAME';
Config.SubP{9,2} = 'D:\enrique_GD_RD\ORT\modelos_cvx\distflow\debug\distflow_whole.log';
Config.SubP{10,1} = 'MSK_IPAR_INFEAS_REPORT_LEVEL';
Config.SubP{10,2} = 20;
Config.SubP{11,1} = 'MSK_IPAR_PRESOLVE_USE';
Config.SubP{11,2} = 'MSK_PRESOLVE_MODE_ON';
Config.SubP{12,1} = 'MSK_IPAR_PRESOLVE_ELIMINATOR_USE';
Config.SubP{12,2} = 'MSK_ON';

% Centralizado
Config.Centr{1, 1} = 'MSK_DPAR_MIO_TOL_REL_RELAX_INT';
Config.Centr{1, 2} = 0.01;
Config.Centr{2, 1} = 'MSK_DPAR_MIO_DISABLE_TERM_TIME';
Config.Centr{2, 2} = TopSecs*TopIT;
Config.Centr{3, 1} = 'MSK_DPAR_MIO_NEAR_TOL_REL_GAP';
Config.Centr{3, 2} = 0.025;
Config.Centr{4, 1} = 'MSK_IPAR_OPF_WRITE_PROBLEM';
Config.Centr{4, 2} = 'MSK_ON';
Config.Centr{5, 1} = 'MSK_IPAR_OPF_WRITE_HEADER';
Config.Centr{5, 2} = 'MSK_ON';
Config.Centr{6, 1} = 'MSK_IPAR_INFEAS_REPORT_AUTO';
Config.Centr{6, 2} = 'MSK_ON';
Config.Centr{7, 1} = 'MSK_SPAR_PARAM_WRITE_FILE_NAME';
Config.Centr{7, 2} = 'D:\enrique_GD_RD\ORT\modelos_cvx\distflow\debug\distflow_whole.log';
Config.Centr{8, 1} = 'MSK_IPAR_INFEAS_REPORT_LEVEL';
Config.Centr{8, 2} = 5;
Config.Centr{9, 1} = 'MSK_IPAR_PRESOLVE_USE';
Config.Centr{9, 2} = 'MSK_PRESOLVE_MODE_ON';
Config.Centr{10, 1} = 'MSK_IPAR_PRESOLVE_ELIMINATOR_USE';
Config.Centr{10, 2} = 'MSK_ON';

% Centralizado inicial
Config.Ini{1, 1} = 'MSK_DPAR_MIO_TOL_REL_RELAX_INT';
Config.Ini{1, 2} = 0.01;
Config.Ini{2, 1} = 'MSK_DPAR_MIO_DISABLE_TERM_TIME';
Config.Ini{2, 2} = TopSecs;
Config.Ini{3, 1} = 'MSK_DPAR_MIO_NEAR_TOL_REL_GAP';
Config.Ini{3, 2} = 0.025;
Config.Ini{4, 1} = 'MSK_IPAR_OPF_WRITE_PROBLEM';
Config.Ini{4, 2} = 'MSK_ON';
Config.Ini{5, 1} = 'MSK_IPAR_OPF_WRITE_HEADER';
Config.Ini{5, 2} = 'MSK_ON';
Config.Ini{6, 1} = 'MSK_IPAR_INFEAS_REPORT_AUTO';
Config.Ini{6, 2} = 'MSK_ON';
Config.Ini{7, 1} = 'MSK_SPAR_PARAM_WRITE_FILE_NAME';
Config.Ini{7, 2} = 'D:\enrique_GD_RD\ORT\modelos_cvx\distflow\debug\distflow_whole.log';
Config.Ini{8, 1} = 'MSK_IPAR_INFEAS_REPORT_LEVEL';
Config.Ini{8, 2} = 5;
Config.Ini{9, 1} = 'MSK_IPAR_PRESOLVE_USE';
Config.Ini{9, 2} = 'MSK_PRESOLVE_MODE_ON';
Config.Ini{10, 1} = 'MSK_IPAR_PRESOLVE_ELIMINATOR_USE';
Config.Ini{10, 2} = 'MSK_ON';
% Subproblemas distribuidos
Config.SubP{1,1} = 'MSK_DPAR_MIO_MAX_TIME';
Config.SubP{1,2} = TopSecs/2;
Config.SubP{2,1} = 'MSK_DPAR_OPTIMIZER_MAX_TIME';
Config.SubP{2,2} = TopSecs;
Config.SubP{3,1} = 'MSK_DPAR_MIO_TOL_REL_RELAX_INT';
Config.SubP{3,2} = 0.01;
Config.SubP{4,1} = 'MSK_DPAR_MIO_DISABLE_TERM_TIME';
Config.SubP{4,2} = 10;
Config.SubP{5,1} = 'MSK_DPAR_MIO_NEAR_TOL_REL_GAP';
Config.SubP{5,2} = 0.025;
Config.SubP{6,1} = 'MSK_IPAR_OPF_WRITE_PROBLEM';
Config.SubP{6,2} = 'MSK_ON';
Config.SubP{7,1} = 'MSK_IPAR_OPF_WRITE_HEADER';
Config.SubP{7,2} = 'MSK_ON';
Config.SubP{8,1} = 'MSK_IPAR_INFEAS_REPORT_AUTO';
Config.SubP{8,2} = 'MSK_ON';
Config.SubP{9,1} = 'MSK_SPAR_PARAM_WRITE_FILE_NAME';
Config.SubP{9,2} = 'D:\enrique_GD_RD\ORT\modelos_cvx\distflow\debug\distflow_whole.log';
Config.SubP{10,1} = 'MSK_IPAR_INFEAS_REPORT_LEVEL';
Config.SubP{10,2} = 20;
Config.SubP{11,1} = 'MSK_IPAR_PRESOLVE_USE';
Config.SubP{11,2} = 'MSK_PRESOLVE_MODE_ON';
Config.SubP{12,1} = 'MSK_IPAR_PRESOLVE_ELIMINATOR_USE';
Config.SubP{12,2} = 'MSK_ON';

% Centralizado
Config.Centr{1, 1} = 'MSK_DPAR_MIO_TOL_REL_RELAX_INT';
Config.Centr{1, 2} = 0.01;
Config.Centr{2, 1} = 'MSK_DPAR_MIO_DISABLE_TERM_TIME';
Config.Centr{2, 2} = TopSecs*TopIT;
Config.Centr{3, 1} = 'MSK_DPAR_MIO_NEAR_TOL_REL_GAP';
Config.Centr{3, 2} = 0.025;
Config.Centr{4, 1} = 'MSK_IPAR_OPF_WRITE_PROBLEM';
Config.Centr{4, 2} = 'MSK_ON';
Config.Centr{5, 1} = 'MSK_IPAR_OPF_WRITE_HEADER';
Config.Centr{5, 2} = 'MSK_ON';
Config.Centr{6, 1} = 'MSK_IPAR_INFEAS_REPORT_AUTO';
Config.Centr{6, 2} = 'MSK_ON';
Config.Centr{7, 1} = 'MSK_SPAR_PARAM_WRITE_FILE_NAME';
Config.Centr{7, 2} = 'D:\enrique_GD_RD\ORT\modelos_cvx\distflow\debug\distflow_whole.log';
Config.Centr{8, 1} = 'MSK_IPAR_INFEAS_REPORT_LEVEL';
Config.Centr{8, 2} = 5;
Config.Centr{9, 1} = 'MSK_IPAR_PRESOLVE_USE';
Config.Centr{9, 2} = 'MSK_PRESOLVE_MODE_ON';
Config.Centr{10, 1} = 'MSK_IPAR_PRESOLVE_ELIMINATOR_USE';
Config.Centr{10, 2} = 'MSK_ON';

% Centralizado inicial
Config.Ini{1, 1} = 'MSK_DPAR_MIO_TOL_REL_RELAX_INT';
Config.Ini{1, 2} = 0.01;
Config.Ini{2, 1} = 'MSK_DPAR_MIO_DISABLE_TERM_TIME';
Config.Ini{2, 2} = TopSecs;
Config.Ini{3, 1} = 'MSK_DPAR_MIO_NEAR_TOL_REL_GAP';
Config.Ini{3, 2} = 0.025;
Config.Ini{4, 1} = 'MSK_IPAR_OPF_WRITE_PROBLEM';
Config.Ini{4, 2} = 'MSK_ON';
Config.Ini{5, 1} = 'MSK_IPAR_OPF_WRITE_HEADER';
Config.Ini{5, 2} = 'MSK_ON';
Config.Ini{6, 1} = 'MSK_IPAR_INFEAS_REPORT_AUTO';
Config.Ini{6, 2} = 'MSK_ON';
Config.Ini{7, 1} = 'MSK_SPAR_PARAM_WRITE_FILE_NAME';
Config.Ini{7, 2} = 'D:\enrique_GD_RD\ORT\modelos_cvx\distflow\debug\distflow_whole.log';
Config.Ini{8, 1} = 'MSK_IPAR_INFEAS_REPORT_LEVEL';
Config.Ini{8, 2} = 5;
Config.Ini{9, 1} = 'MSK_IPAR_PRESOLVE_USE';
Config.Ini{9, 2} = 'MSK_PRESOLVE_MODE_ON';
Config.Ini{10, 1} = 'MSK_IPAR_PRESOLVE_ELIMINATOR_USE';
Config.Ini{10, 2} = 'MSK_ON';


%% Llamada al modelo

cantTaps = length(Trafos);
cantCaps = length(Caps);
cantCargs = length(Cargas);

Data.Red.Bus.uTop = Data.Red.Bus.uTop*1.5;
Data.Red.Bus.uLow = Data.Red.Bus.uLow/1.5;
Data.Red.Bus.uTop(1,:) = 1;
Data.Red.Bus.uLow(1,:) = 1;

[Var_nxn, opt_nxn, status, Dmod] = llamarCentralizado(Data, Config, utilCarg);

Data.Red.Bus.pCLow = full(Data.Red.Bus.pCLow);
Data.Red.Bus.qCLow = full(Data.Red.Bus.qCLow);

Data.Gen.Tras.pgLow = full(Data.Gen.Tras.pgLow);
Data.Gen.Tras.qgLow = full(Data.Gen.Tras.qgLow);
Data.Gen.Tras.pgTop = full(Data.Gen.Tras.pgTop);
Data.Gen.Tras.qgTop = full(Data.Gen.Tras.qgTop);

Data.Cost.cdv = repmat(full(Data.Cost.cdv), [1 Config.Etapas]);

Data.Red.Branch.r = repmat(full(Data.Red.Branch.r), [1 1 Config.Etapas]);
Data.Red.Branch.x = repmat(full(Data.Red.Branch.x), [1 1 Config.Etapas]);
Data.Red.Branch.lTop = repmat(full(Data.Red.Branch.lTop), [1 1 Config.Etapas]);
Data.Red.Branch.yTop = repmat(full(Data.Red.Branch.yTop), [1 1 Config.Etapas]);
Data.Red.Branch.yLow = repmat(full(Data.Red.Branch.yLow), [1 1 Config.Etapas]);


Data.Red.Bus.Ntr = repmat(Data.Red.Bus.Ntr, [1 Config.Etapas]);
Data.Red.Bus.uLow = repmat(Data.Red.Bus.uLow, [1 Config.Etapas]);
Data.Red.Bus.uTop = repmat(Data.Red.Bus.uTop, [1 Config.Etapas]);
Data.Red.Bus.TapLow = repmat(Data.Red.Bus.TapLow, [1 Config.Etapas]);
Data.Red.Bus.TapTop = repmat(Data.Red.Bus.TapTop, [1 Config.Etapas]);
Data.Red.Bus.indTap = repmat(Data.Red.Bus.indTap, [1 Config.Etapas]);

Data.Red.Bus.Ncp = repmat(Data.Red.Bus.Ncp, [1 Config.Etapas]);
Data.Red.Bus.CapLow = repmat(Data.Red.Bus.CapLow, [1 Config.Etapas]);
Data.Red.Bus.CapTop = repmat(Data.Red.Bus.CapTop, [1 Config.Etapas]);
Data.Red.Bus.indCap = repmat(Data.Red.Bus.indCap, [1 Config.Etapas]);

Data.Util.pzCnLowE = repmat(Data.Util.pzCnLowE, [1 Config.Etapas]);
Data.Util.pzCnTopE = repmat(Data.Util.pzCnTopE, [1 Config.Etapas]);
Data.Util.qzCnLowE = repmat(Data.Util.qzCnLowE, [1 Config.Etapas]);
Data.Util.qzCnTopE = repmat(Data.Util.qzCnTopE, [1 Config.Etapas]);
Data.Util.pzCnPrefE = repmat(Data.Util.pzCnPrefE, [1 Config.Etapas]);

Data.Util.pzCnPref = Data.Util.pzCnPref(:,(1:Config.Etapas),:);
Data.Util.pzCnLow = Data.Util.pzCnLow(:,(1:Config.Etapas),:);
Data.Util.pzCnTop = Data.Util.pzCnTop(:,(1:Config.Etapas),:);

% Data.Util.aE = repmat(Data.Util.aE, [1 Config.Etapas]);

Data.Util.betaT = repmat(Data.Util.betaT, [1 Config.Etapas]);

a1 = repmat(Data.Red.Bus.alpha(:,1), [1 Config.Etapas 2]);
a1(:,:,2) = repmat(Data.Red.Bus.alpha(:,2), [1 Config.Etapas]);
Data.Red.Bus.alpha = a1;

Data.Red.Bus.pCLowm = Data.Red.Bus.pCLow(:,(1:Config.Etapas));
Data.Red.Bus.qCLowm = Data.Red.Bus.qCLow(:,(1:Config.Etapas));

Data.Gen.Tras.pgLowm = Data.Gen.Tras.pgLow(:,(1:Config.Etapas));
Data.Gen.Tras.qgLowm = Data.Gen.Tras.qgLow(:,(1:Config.Etapas));
Data.Gen.Tras.pgTopm = Data.Gen.Tras.pgTop(:,(1:Config.Etapas));
Data.Gen.Tras.qgTopm = Data.Gen.Tras.qgTop(:,(1:Config.Etapas));

Data.Cost.m = Data.Cost.m(:,(1:Config.Etapas));

Data.Cost.piPTrasm = Data.Cost.piPTras(:,(1:Config.Etapas));
Data.Cost.piQmtras = Data.Cost.piQmtras(:,(1:Config.Etapas));
Data.Cost.piQMtras = Data.Cost.piQMtras(:,(1:Config.Etapas));

Data.Cost.cdvm = Data.Cost.cdv;


Data.Util.betaE = Data.Util.betaE(:,(1:Config.Etapas));



Data.St.AC.epsilon = repmat(Data.St.AC.epsilon, [1 Config.Etapas]);
Data.temp = repmat(Data.temp, [1 Config.Etapas]);
Data.St.AC.eta = repmat(Data.St.AC.eta, [1 Config.Etapas]);
Data.St.AC.tempLow = repmat(Data.St.AC.tempLow, [1 Config.Etapas]);
Data.St.AC.tempTop = repmat(Data.St.AC.tempTop, [1 Config.Etapas]);



[Var_m, opt_m] = distflowCentralizadoM(Data, Config);

printSalidasDistflow(Var_nxn, Dmod, Config, cantTaps, cantCaps, cantCargs, outFilename_c, [], [], [], [], []);
printSalidasDistflow(Var_m, Dmod, Config, cantTaps, cantCaps, cantCargs, outFilename_r, [], [], [], [], []);

