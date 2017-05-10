function [Data] = Dfig_500kw(Data, windNodes)
%% Datos Generador Eolico
w = 1; % Omega? o ws? TODO sacar
ws = pi*100;

% Np?
Np = 2;
% G?
G = 90;
% R_ = 1;
R_ = 37.5;
% lambda_opt ?
lambda_opt = 6.325;
Omega = 169.6464;

% rho = 1;
% N_er = 1;
rho = 1.225;
N_er = 1;

% % % RW_12 = 0.00000052939453;
% % % LW_12 = 0.00000031992188;
% % RW_12 = 3.66211e-5;
% % LW_12 = .019174805;
RW_12 = .014648438;
LW_12 = .030492541;
% RW_12 = .01375;
% LW_12 = .23561925;

% RW_13 = 0.00000052939453;
% LW_13 = 0.00000031992188;
RW_13 = 0.03;
LW_13 = 0.3;

% % % RW_45 = 0.00000052939453;
% % % LW_45 = 0.00000031992188;
% % RW_45 = 4.88281e-5;
% % LW_45 = .019174805;
RW_45 = .017407227;
LW_45 = .030492541;
% RW_45 = .011875;
% LW_45 = .09817469;

lW_12_Top = 3;
lW_13_Top = 3;
lW_45_Top = 3;
sW_3_Top = .2;
sW_5_Top = .6;
% sW_5_Top = .2;

n_ = .95;

vmm = 5;
% vm = 11.1755922;
vm = 9.55;
vM = 11.90354704;
vMM = 25;

cvW_3 = .01;
cvW_5 = .01;
crW_3 = .03;
crW_5 = .03;

P_nMec =5e5;
c1=0.055;
c2=116;
c3=0.4;
c4=5;
c5=12.5;
c6=0.08;
c7=0.035;

uLow_2 = sqrt(.49);
uTop_2 = 1.1;

uLow_3 = sqrt(.49);
uTop_3 = 1.1;

PQnorm_2 = 2;
PQnorm_3 = .5;

% TODO sacar
C_plb = 1;
vv = 1;

% windNodes = [4; 7; 3];
% windNodes = [4];
% windNodes = [];
lenWN = length(windNodes);
    
% P_mec = rho*pi*R_^2*C_plb*vv^3;
% P_mec = 0.15271;

Data.Gen.DFIG.Tg = zeros(5,5);
Data.Gen.DFIG.Tg(1,2) = 1;
Data.Gen.DFIG.Tg(1,3) = 1;
Data.Gen.DFIG.Tg(4,5) = 1;

% Data.Gen.Wa.Sg = [3];
% indSg = zeros(1,5);
% indSg(Data.Gen.Wa.Sg) = 1;

Data.Gen.DFIG.r = zeros(5,5,lenWN);
Data.Gen.DFIG.x = Data.Gen.DFIG.r;
Data.Gen.DFIG.lTop = Data.Gen.DFIG.r;

Data.Gen.DFIG.cv = zeros(5,1,lenWN);
Data.Gen.DFIG.cr = Data.Gen.DFIG.cv;
Data.Gen.DFIG.sTop = Data.Gen.DFIG.cv;
Data.Gen.DFIG.xiTop = Data.Gen.DFIG.cv;
Data.Gen.DFIG.I = zeros(size(Data.Gen.Tras.I));
Data.Gen.DFIG.I(windNodes) = 1;

Data.Gen.DFIG.n_ = n_;
Data.Gen.DFIG.N_er = N_er;
Data.Gen.DFIG.w = w;
Data.Gen.DFIG.rho = rho;
Data.Gen.DFIG.R_ = R_;
Data.Gen.DFIG.C_plb = C_plb;
Data.Gen.DFIG.vv = vv;
Data.Gen.DFIG.P_mec = 0;

Data.Gen.DFIG.vmm = vmm;
Data.Gen.DFIG.vm = vm;
Data.Gen.DFIG.vM = vM;
Data.Gen.DFIG.vMM = vMM;

Data.Gen.DFIG.P_nMec = P_nMec;
Data.Gen.DFIG.c1 = c1;
Data.Gen.DFIG.c2 = c2;
Data.Gen.DFIG.c3 = c3;
Data.Gen.DFIG.c4 = c4;
Data.Gen.DFIG.c5 = c5;
Data.Gen.DFIG.c6 = c6;
Data.Gen.DFIG.c7 = c7;

Data.Gen.DFIG.Omega = Omega;
Data.Gen.DFIG.G = G;
Data.Gen.DFIG.Np = Np;
Data.Gen.DFIG.rho = rho;
Data.Gen.DFIG.ws = ws;
Data.Gen.DFIG.lambda_opt = lambda_opt;

Data.Gen.DFIG.uLow = Data.Gen.DFIG.cr;
Data.Gen.DFIG.uTop = ones(size(Data.Gen.DFIG.uLow))*3;
% Data.Gen.DFIG.PQnorm = zeros(size(Data.Gen.DFIG.lTop));

% for wnd = 1:lenWN
	Data.Gen.DFIG.rIE = RW_12;
    Data.Gen.DFIG.rIF = RW_13;
    Data.Gen.DFIG.rOR = RW_45;

	Data.Gen.DFIG.xIE = Data.Gen.DFIG.w*LW_12;
    Data.Gen.DFIG.xIF = Data.Gen.DFIG.w*LW_13;
    Data.Gen.DFIG.xOR = Data.Gen.DFIG.w*LW_45;

	Data.Gen.DFIG.lTopIE = lW_12_Top;
    Data.Gen.DFIG.lTopIF = lW_13_Top;
    Data.Gen.DFIG.lTopOR = lW_45_Top;

	Data.Gen.DFIG.cvF = cvW_3;
    Data.Gen.DFIG.cvR = cvW_5;

	Data.Gen.DFIG.crF = crW_3;
    Data.Gen.DFIG.crR = crW_5;

	Data.Gen.DFIG.uLowE = uLow_2;
    Data.Gen.DFIG.uLowF = uLow_3;

    Data.Gen.DFIG.uTopE = uTop_2;
	Data.Gen.DFIG.uTopF = uTop_3;

	Data.Gen.DFIG.sTopF = sW_3_Top;
    Data.Gen.DFIG.sTopR = sW_5_Top;

	Data.Gen.DFIG.xiTopF = sW_3_Top^2;
    Data.Gen.DFIG.xiTopR = sW_5_Top^2;

	Data.Gen.DFIG.Ig(3,1,:) = 1;

    Data.Gen.DFIG.PQnormIE = PQnorm_2;
    Data.Gen.DFIG.PQnormIF = PQnorm_3;
    
    
%     ini = 5*(wnd-1)+1;
%     fin = 5*(wnd-1)+5;
% 
%     Data.Gen.DFIG.r(:,ini:fin) = Data.Gen.DFIG.r(:,ini:fin) .* Data.Gen.DFIG.Tg;
% 
%     Data.Gen.DFIG.x(:,ini:fin) = Data.Gen.DFIG.x(:,ini:fin) .* Data.Gen.DFIG.Tg;
% 
%     Data.Gen.DFIG.lTop(:,ini:fin) = Data.Gen.DFIG.lTop(:,ini:fin) .* Data.Gen.DFIG.Tg;
% 
%     Data.Gen.DFIG.cv(ini:fin) = Data.Gen.DFIG.cv(ini:fin) .* indSg;
% 
% 	Data.Gen.DFIG.cr(ini:fin) = Data.Gen.DFIG.cr(ini:fin) .* indSg;
% 
% 	Data.Gen.DFIG.sTop(ini:fin) = Data.Gen.DFIG.sTop(ini:fin) .* indSg;
% 
% 	Data.Gen.DFIG.Ig(ini:fin) = Data.Gen.DFIG.Ig(ini:fin) .* indSg;

    
% end
