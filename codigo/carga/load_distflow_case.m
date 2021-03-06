function [Data] = load_distflow_case(Data, xls_file, bus_sheet, branch_sheet, Trafos, Caps, Cargas, App, Switches)



[n_br,~,~] = xlsread(xls_file, branch_sheet);
[n_bu,~,~] = xlsread(xls_file, bus_sheet);

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

Data.Red.Branch.T(:,:) = T+T';
Data.Red.Branch.r(:,:) = r+r';
Data.Red.Branch.x(:,:) = x+x';
Data.Red.Branch.lTop(:,:) = lTop+lTop';
Data.Red.Branch.yTop(:,:) = Data.Red.Branch.T;
Data.Red.Branch.yLow(:,:) = Data.Red.Branch.T;
Data.Red.Branch.Tswitches(:,:) = Data.Red.Branch.T*0;

if Switches.all
    % Todos los arcos son decidibles
	Data.Red.Branch.Tswitches(:,:) = Data.Red.Branch.T;
    Data.Red.Branch.yLow(:,:) = Data.Red.Branch.T*0;
end
    
for i=1:length(Switches.i)
	Data.Red.Branch.yLow(Switches.i(i), Switches.j(i))=0;
	Data.Red.Branch.yLow(Switches.j(i), Switches.i(i))=0;

	Data.Red.Branch.Tswitches(Switches.i(i), Switches.j(i))=1;
	Data.Red.Branch.Tswitches(Switches.j(i), Switches.i(i))=1;
end

Data.Red.Bus.uLow(:,:) = sparse(n_bu(:,5));
Data.Red.Bus.uTop(:,:) = sparse(n_bu(:,4));

for a = 1:length(App)
    Data.Red.Bus.alpha(:,App(a).I) = App(a).alpha;
end


%% Transformadores
for estT = 1:size(Trafos,1)
	Data.Red.Branch.NtrLow(Trafos(estT).nodI, Trafos(estT).nodJ) = min(Trafos(estT).N);
	Data.Red.Branch.NtrTop(Trafos(estT).nodI, Trafos(estT).nodJ) = max(Trafos(estT).N);
	Data.Red.Branch.NtrIni(Trafos(estT).nodI, Trafos(estT).nodJ) = Trafos(estT).ini;
    if Trafos(estT).reg
        Data.Red.Branch.Itreg(Trafos(estT).nodI, Trafos(estT).nodJ) = 1;
    else
        Data.Red.Branch.Itap(Trafos(estT).nodI, Trafos(estT).nodJ) = 1;
    end
	Data.Red.Branch.Tap(Trafos(estT).nodI, Trafos(estT).nodJ) = Trafos(estT).TP;
end

Data.Red.Branch.NtrLow = Data.Red.Branch.NtrLow + Data.Red.Branch.NtrLow';
Data.Red.Branch.NtrTop = Data.Red.Branch.NtrTop + Data.Red.Branch.NtrTop';
Data.Red.Branch.NtrIni = Data.Red.Branch.NtrIni + Data.Red.Branch.NtrIni';
Data.Red.Branch.Itap = Data.Red.Branch.Itap + Data.Red.Branch.Itap';
Data.Red.Branch.Tap  = Data.Red.Branch.Tap + Data.Red.Branch.Tap';


%% Capacitores
for estC = 1:size(Caps,1)
	Data.Red.Bus.NcpLow(Caps(estC).nod) = min(Caps(estC).N);
	Data.Red.Bus.NcpTop(Caps(estC).nod) = max(Caps(estC).N);
	Data.Red.Bus.NcpIni(Caps(estC).nod) = Caps(estC).ini;
	Data.Red.Bus.Icap(Caps(estC).nod) = 1;
	Data.Red.Bus.Cap(Caps(estC).nod) = Caps(estC).TP;
end

%% Cargas no interrumpibles
for estCg = 1:length(Cargas)
	Data.ClNI.pC(Cargas(estCg).nod) = Cargas(estCg).pC;
	Data.ClNI.qC(Cargas(estCg).nod) = Cargas(estCg).qC;
	Data.ClNI.d(Cargas(estCg).nod) = Cargas(estCg).dur;
	Data.ClNI.I(Cargas(estCg).nod) = 1;
    Data.ClNI.nMultipTop(Cargas(estCg).nod) = Cargas(estCg).nMultipTop;
    Data.ClNI.nMultipLow(Cargas(estCg).nod) = Cargas(estCg).nMultipLow;
end
