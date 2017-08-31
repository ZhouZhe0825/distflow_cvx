function [Data] = cargaCostosDefault(Data, Trafos, Caps, Switches, GenBas, fileCostosTension, Tras, Solares, Eolicos)

[Data.Cost.m, Data.Cost.delta, Data.Cost.cdv] = costosTension(Data,fileCostosTension);

[Data.Cost.piPTras, Data.Cost.piQmtras, Data.Cost.piQMtras] = costosTrasmision(Data,Tras);

[Data.Cost.rhopPv, Data.Cost.rhomqPv, Data.Cost.rhoMqPv] = costosPv(Data,Solares);

[Data.Cost.rhopWi, Data.Cost.rhomqWi, Data.Cost.rhoMqWi] = costosDfig(Data,Eolicos);

for i = 1:length(Caps)
	Data.Cost.cCap(Caps(i).nod) = Caps(i).cambio;
end    

for i = 1:length(Trafos)
	Data.Cost.cTap(Trafos(i).nodI,Trafos(i).nodJ) = Trafos(i).cambio;
end    

Data.Cost.cTap = Data.Cost.cTap + Data.Cost.cTap';

Data.Cost.cY = Data.Red.Branch.T .* Switches.cY;

[Data.Cost.cBas] = costosGBasico(Data,GenBas);
