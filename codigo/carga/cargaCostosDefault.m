function [Data] = cargaCostosDefault(Data)

Data.Cost.rhopPv = ones(size(Data.Gen.Pv.I)).*Data.Gen.Pv.I;
Data.Cost.rhomqPv = ones(size(Data.Gen.Pv.I)).*Data.Gen.Pv.I;
Data.Cost.rhoMqPv = ones(size(Data.Gen.Pv.I)).*Data.Gen.Pv.I;

Data.Cost.piPTras = zeros(size(Data.Gen.Pv.I));
Data.Cost.piQmtras = zeros(size(Data.Gen.Pv.I));
Data.Cost.piQMtras = zeros(size(Data.Gen.Pv.I));

Data.Cost.piPTras(Data.Red.Bus.v0) = 1;
Data.Cost.piQmtras(Data.Red.Bus.v0) = 1;
Data.Cost.piQMtras(Data.Red.Bus.v0) = 1; 

Data.Cost.rhopWi = zeros(length(Data.Red.Branch.T),1);
Data.Cost.rhomqWi = zeros(length(Data.Red.Branch.T),1);
Data.Cost.rhoMqWi = zeros(length(Data.Red.Branch.T),1);
Data.Cost.cdv = ones(length(Data.Red.Branch.T),1);
