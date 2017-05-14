function [Data] = cargaCostosDefault(Data, Trafos, Caps, Switches, mHor, cdv, delta, rhopPv, rhomqPv, rhoMqPv, rhopWi, rhomqWi, rhoMqWi, piPTrasHor, piQmtrasHor, piQMtrasHor)

Data.Cost.rhopPv = rhopPv;
Data.Cost.rhomqPv = rhomqPv;
Data.Cost.rhoMqPv = rhoMqPv;

Data.Cost.rhopWi = rhopWi;
Data.Cost.rhomqWi = rhomqWi;
Data.Cost.rhoMqWi = rhoMqWi;

Data.Cost.piPTras = Data.Gen.Tras.I * piPTrasHor';
Data.Cost.piQmtras = Data.Gen.Tras.I * piQmtrasHor';
Data.Cost.piQMtras = Data.Gen.Tras.I * piQMtrasHor';

Data.Cost.m = mHor';
Data.Cost.delta = delta;
Data.Cost.cdv = cdv;

Data.Cost.cCap = Data.Red.Bus.indCap;
for i = 1:length(Caps)
	Data.Cost.cCap(Caps(i).nod) = Caps(i).cambio;
end    

Data.Cost.cTap = Data.Red.Bus.indTap;
for i = 1:length(Trafos)
	Data.Cost.cTap(Trafos(i).nod) = Trafos(i).cambio;
end    

Data.Cost.cY = Data.Red.Branch.T .* Switches.cY;
