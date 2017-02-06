function [Data] = load_distflow_case(xls_file, bus_sheet, branch_sheet)



[n_br,t_br,r_br] = xlsread(xls_file, branch_sheet);
[n_bu,t_bu,r_bu] = xlsread(xls_file, bus_sheet);

% Por ahora funciona solo para red pasiva
n = length(n_bu(:,1));
i = n_br(:,1);
j = n_br(:,2);
T_br = ones(length(n_br(:,1)),1);
R_br = n_br(:,3);
X_br = n_br(:,4);
lTop_br = n_br(:,5);
T = sparse(i, j, T_br, n, n);
r = sparse(i, j, R_br, n, n);
x = sparse(i, j, X_br, n, n);
lTop = sparse(i, j, lTop_br, n, n);
V_bu = ones(n,1);

Data.Red.Branch.T = T+T';
Data.Red.Branch.r = r+r';
Data.Red.Branch.x = x+x';
Data.Red.Branch.lTop = lTop+lTop';

Data.Red.Branch.cY = Data.Red.Branch.T;

Data.Red.Bus.alpha = sparse(1 * V_bu);
Data.Red.Bus.uLow = sparse(n_bu(:,7));
Data.Red.Bus.uTop = sparse(n_bu(:,6));
Data.Red.Bus.pCLow = sparse(n_bu(:,3));
Data.Red.Bus.qCLow = sparse(n_bu(:,4));
Data.Red.Bus.Ncp = zeros(n, 1);
Data.Red.Bus.v0 = find(~isnan(n_bu(:,2)));
Data.Red.Bus.Q0Top = 10;
Data.Red.Bus.Q0Low = -10;
Data.Red.Bus.P0Top = 10;
Data.Red.Bus.P0Low = -10;
Data.Red.Bus.Ntr = zeros(n, 1);

% Data.Red.Cost = 0;

Data.Gen.Pv.I = sparse(0 * V_bu);
Data.Gen.Pv.cv = sparse(0 * V_bu);
Data.Gen.Pv.cr = sparse(0 * V_bu);
Data.Gen.Pv.qgTop = sparse(0 * V_bu);
Data.Gen.Pv.sTop = sparse(0 * V_bu);
Data.Gen.Pv.pgTop = sparse(0 * V_bu);


% T = T+T';
% r = r+r';
% r = x+x';
% lTop = lTop+lTop';
% V_bu = ones(n,1);
% I = 0 * V_bu;
% alpha = 0.4 * V_bu;
% cv = 0 * V_bu;
% cr = 0 * V_bu;
% qgTop = 0 * V_bu;
% uLow = n_bu(:,7);
% uTop = n_bu(:,6);
% pcLow = n_bu(:,3);
% qcLow = n_bu(:,4);
% qcCap = 0 * V_bu;
% sTop = 0 * V_bu;
% pgTop = 0 * V_bu;
% v0 = find(~isnan(n_bu(:,2)));
% Q0Top = 10;
% Q0Low = -10;
% P0Top = 10;
% P0Low = -10;
