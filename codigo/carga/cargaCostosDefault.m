function [Data] = cargaCostosDefault(Data, Trafos, Caps, Switches, fileCostosTension, fileCostosTras, Solares, Eolicos)

[Data.Cost.m, Data.Cost.delta, Data.Cost.cdv] = costosTension(Data,fileCostosTension);

[Data.Cost.piPTras, Data.Cost.piQmtras, Data.Cost.piQMtras] = costosTrasmision(Data,fileCostosTras);

[Data.Cost.rhopPv, Data.Cost.rhomqPv, Data.Cost.rhoMqPv] = costosPv(Data,Solares);

[Data.Cost.rhopWi, Data.Cost.rhomqWi, Data.Cost.rhoMqWi] = costosDfig(Data,Eolicos);

Data.Cost.cCap = Data.Red.Bus.Icap;
for i = 1:length(Caps)
	Data.Cost.cCap(Caps(i).nod) = Caps(i).cambio;
end    

for i = 1:length(Trafos)
	Data.Cost.cTap(Trafos(i).nodI,Trafos(i).nodJ) = Trafos(i).cambio;
end    

Data.Cost.cTap = Data.Cost.cTap + Data.Cost.cTap';

Data.Cost.cY = Data.Red.Branch.T .* Switches.cY;
