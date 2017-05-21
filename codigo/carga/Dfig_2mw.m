function [w,ws,Np,G,R_,lambda_opt,Omega,rho,N_er,RW_12,LW_12,RW_13,LW_13,RW_45,LW_45,lW_12_Top,lW_13_Top,lW_45_Top,sW_3_Top,sW_5_Top,vmm,vm,vM,vMM,cvW_3,cvW_5,crW_3,crW_5,P_nMec,c1,c2,c3,c4,c5,c6,c7,uLow_2,uTop_2,uLow_3,uTop_3,PQnorm_2,PQnorm_3,C_plb,vv] = Dfig_2mw()
%% Datos Generador Eolico
w = 1; % Omega? o ws? TODO sacar
ws = pi*100;

Np = 2;
G = 90;
R_ = 37.5;
lambda_opt = 6.325;
Omega = 169.6464;

rho = 1.225;
N_er = 1;

RW_12 = .014648438;
LW_12 = .030492541;

RW_13 = 0.03;
LW_13 = 0.3;

RW_45 = .017407227;
LW_45 = .030492541;

lW_12_Top = 3;
lW_13_Top = 3;
lW_45_Top = 3;
sW_3_Top = 3;
sW_5_Top = 3;

vmm = 5;
vm = 9.55;
vM = 11.90354704;
vMM = 25;

cvW_3 = .01;
cvW_5 = .01;
crW_3 = .03;
crW_5 = .03;

P_nMec =2e6;
c1=0.22;
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

C_plb = 1;
vv = 1;

