function [Data] = Dfig_250kw(Data, windNodes)
%% Datos Generador Eolico
w = 1; % Omega? o ws? TODO sacar
ws = pi*100;

et = size(Data.Gen.DFIG.rIE,2);

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
sW_3_Top = 3;
sW_5_Top = 3;
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

P_nMec =25e5;
c1=0.0275;
c2=116;
c3=0.4;
c4=5;
c5=12.5;
c6=0.08;
c7=0.035;

uLow_2 = sqrt(.49);
uTop_2 = sqrt(1.21550625);

uLow_3 = sqrt(.49);
uTop_3 = sqrt(1.21550625);

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


Data.Gen.DFIG.n_(windNodes) = n_;
Data.Gen.DFIG.N_er(windNodes) = N_er;
Data.Gen.DFIG.w(windNodes) = w;
Data.Gen.DFIG.rho(windNodes) = rho;
Data.Gen.DFIG.R_(windNodes) = R_;
Data.Gen.DFIG.C_plb(windNodes) = C_plb;
Data.Gen.DFIG.vv(windNodes) = vv;
Data.Gen.DFIG.P_mec(windNodes) = 0;

Data.Gen.DFIG.vmm(windNodes) = vmm;
Data.Gen.DFIG.vm(windNodes) = vm;
Data.Gen.DFIG.vM(windNodes) = vM;
Data.Gen.DFIG.vMM(windNodes) = vMM;

Data.Gen.DFIG.P_nMec(windNodes) = P_nMec;
Data.Gen.DFIG.c1(windNodes) = c1;
Data.Gen.DFIG.c2(windNodes) = c2;
Data.Gen.DFIG.c3(windNodes) = c3;
Data.Gen.DFIG.c4(windNodes) = c4;
Data.Gen.DFIG.c5(windNodes) = c5;
Data.Gen.DFIG.c6(windNodes) = c6;
Data.Gen.DFIG.c7(windNodes) = c7;

Data.Gen.DFIG.Omega(windNodes) = Omega;
Data.Gen.DFIG.G(windNodes) = G;
Data.Gen.DFIG.Np(windNodes) = Np;
Data.Gen.DFIG.rho(windNodes) = rho;
Data.Gen.DFIG.ws(windNodes) = ws;
Data.Gen.DFIG.lambda_opt(windNodes) = lambda_opt;

% for wnd = 1:lenWN
	Data.Gen.DFIG.rIE(windNodes,:) = RW_12;
    Data.Gen.DFIG.rIF(windNodes,:) = RW_13;
    Data.Gen.DFIG.rOR(windNodes,:) = RW_45;

	Data.Gen.DFIG.xIE(windNodes,:) = repmat(Data.Gen.DFIG.w(windNodes)*LW_12, [1,et]);
    Data.Gen.DFIG.xIF(windNodes,:) = repmat(Data.Gen.DFIG.w(windNodes)*LW_13, [1,et]);
    Data.Gen.DFIG.xOR(windNodes,:) = repmat(Data.Gen.DFIG.w(windNodes)*LW_45, [1,et]);

	Data.Gen.DFIG.lTopIE(windNodes,:) = lW_12_Top;
    Data.Gen.DFIG.lTopIF(windNodes,:) = lW_13_Top;
    Data.Gen.DFIG.lTopOR(windNodes,:) = lW_45_Top;

	Data.Gen.DFIG.cvF(windNodes,:) = cvW_3;
    Data.Gen.DFIG.cvR(windNodes,:) = cvW_5;

	Data.Gen.DFIG.crF(windNodes,:) = crW_3;
    Data.Gen.DFIG.crR(windNodes,:) = crW_5;

	Data.Gen.DFIG.uLowE(windNodes,:) = uLow_2;
    Data.Gen.DFIG.uLowF(windNodes,:) = uLow_3;

    Data.Gen.DFIG.uTopE(windNodes,:) = uTop_2;
	Data.Gen.DFIG.uTopF(windNodes,:) = uTop_3;

	Data.Gen.DFIG.sTopF(windNodes,:) = sW_3_Top;
    Data.Gen.DFIG.sTopR(windNodes,:) = sW_5_Top;

	Data.Gen.DFIG.xiTopF(windNodes,:) = sW_3_Top^2;
    Data.Gen.DFIG.xiTopR(windNodes,:) = sW_5_Top^2;

% 	Data.Gen.DFIG.Ig(3,1,:) = 1;

    Data.Gen.DFIG.I(windNodes,:) = 1;
    
    Data.Gen.DFIG.PQnormIE(windNodes,:) = PQnorm_2;
    Data.Gen.DFIG.PQnormIF(windNodes,:) = PQnorm_3;
    
    
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
