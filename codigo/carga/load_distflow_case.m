function [Data] = load_distflow_case(xls_file, bus_sheet, branch_sheet, Trafos, Caps, Cargas, App, Switches, cambioTap, cambioCap)



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
Data.Red.Branch.yTop = Data.Red.Branch.T;

if Switches.all
    Data.Red.Branch.yLow = Data.Red.Branch.T;
else
    Data.Red.Branch.yLow = Data.Red.Branch.T * 0;
end
    
for i=1:length(Switches.i)
	Data.Red.Branch.yLow(Switches.i(i), Switches.j(i))=0;
	Data.Red.Branch.yLow(Switches.j(i), Switches.i(i))=0;
end

Data.Red.Bus.alpha = sparse(1 * V_bu);
Data.Red.Bus.uLow = sparse(n_bu(:,7));
Data.Red.Bus.uTop = sparse(n_bu(:,6));
Data.Red.Bus.pCLow = sparse(n_bu(:,3));
Data.Red.Bus.qCLow = sparse(n_bu(:,4));
v0 = find(~isnan(n_bu(:,2)));
Q0Top = 10;
Q0Low = -10;
P0Top = 10;
P0Low = -10;
Data.Red.Bus.alpha = repmat(full(Data.Red.Bus.alpha), [1,length(App)]);
for a = 1:length(App)
    Data.Red.Bus.alpha(:,App(a).I) = App(a).alpha;
end


%% Transformadores
Data.Red.Bus.NtrLow = Data.Red.Bus.uLow * 0;
Data.Red.Bus.NtrTop = Data.Red.Bus.NtrLow;
Data.Red.Bus.NtrIni = Data.Red.Bus.NtrLow;
Data.Red.Bus.indTap = Data.Red.Bus.NtrLow;
Data.Red.Bus.Tap = zeros(n, 1);

for estT = 1:size(Trafos,1)
	Data.Red.Bus.NtrLow(Trafos(estT).nod) = min(Trafos(estT).N);
	Data.Red.Bus.NtrTop(Trafos(estT).nod) = max(Trafos(estT).N);
	Data.Red.Bus.NtrIni(Trafos(estT).nod) = Trafos(estT).ini;
	Data.Red.Bus.indTap(Trafos(estT).nod) = 1;
	Data.Red.Bus.Tap(Trafos(estT).nod) = Trafos(estT).TP;
end

%% Capacitores
Data.Red.Bus.NcpLow = Data.Red.Bus.uLow * 0;
Data.Red.Bus.NcpTop = Data.Red.Bus.NcpLow;
Data.Red.Bus.NcpIni = Data.Red.Bus.NcpLow;
Data.Red.Bus.indCap = Data.Red.Bus.NcpLow;
Data.Red.Bus.Cap = zeros(n, 1);


for estC = 1:size(Caps,1)
	Data.Red.Bus.NcpLow(Caps(estC).nod) = min(Caps(estC).N);
	Data.Red.Bus.NcpTop(Caps(estC).nod) = max(Caps(estC).N);
	Data.Red.Bus.NcpIni(Caps(estC).nod) = Caps(estC).ini;
	Data.Red.Bus.indCap(Caps(estC).nod) = 1;
	Data.Red.Bus.Cap(Caps(estC).nod) = Caps(estC).TP;
end

%% Trasmision
Data.Gen.Tras.pgLow = zeros(n, 1);
Data.Gen.Tras.qgLow = Data.Gen.Tras.pgLow;
Data.Gen.Tras.pgTop = Data.Gen.Tras.pgLow;
Data.Gen.Tras.qgTop = Data.Gen.Tras.pgLow;

Data.Gen.Tras.pgLow(v0,:) = P0Low;
Data.Gen.Tras.qgLow(v0,:) = Q0Low;
Data.Gen.Tras.pgTop(v0,:) = P0Top;
Data.Gen.Tras.qgTop(v0,:) = Q0Top;

Data.Gen.Tras.I = zeros(n, 1);
Data.Gen.Tras.I(v0) = 1;


%% Cargas no interrumpibles
Data.ClNI.pC = Data.Red.Bus.uLow * 0;
Data.ClNI.qC = Data.ClNI.pC;
Data.ClNI.d = Data.ClNI.pC;
Data.ClNI.I = Data.ClNI.pC;
Data.ClNI.nMultipTop = Data.ClNI.pC;
Data.ClNI.nMultipLow = Data.ClNI.pC;

for estCg = 1:length(Cargas)
	Data.ClNI.pC(Cargas(estCg).nod) = Cargas(estCg).pC;
	Data.ClNI.qC(Cargas(estCg).nod) = Cargas(estCg).qC;
	Data.ClNI.d(Cargas(estCg).nod) = Cargas(estCg).dur;
	Data.ClNI.I(Cargas(estCg).nod) = 1;
    Data.ClNI.nMultipTop(Cargas(estCg).nod) = Cargas(estCg).nMultipTop;
    Data.ClNI.nMultipLow(Cargas(estCg).nod) = Cargas(estCg).nMultipLow;
end
