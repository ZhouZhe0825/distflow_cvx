function [Data] = cargaCostosDefault(Data, mHor, piPTrasHor, rhoP_ct, rhoQ_ct, m, delta)

Data.Cost.rhopPv = ones(size(Data.Gen.Pv.I)).*Data.Gen.Pv.I;
Data.Cost.rhomqPv = ones(size(Data.Gen.Pv.I)).*Data.Gen.Pv.I;
Data.Cost.rhoMqPv = ones(size(Data.Gen.Pv.I)).*Data.Gen.Pv.I;

Data.Cost.piPTras = zeros(size(Data.Gen.Pv.I));
Data.Cost.piQmtras = zeros(size(Data.Gen.Pv.I));
Data.Cost.piQMtras = zeros(size(Data.Gen.Pv.I));

v0 = find(Data.Gen.Tras.I == 1);
Data.Cost.piPTras(v0) = 1;
Data.Cost.piQmtras(v0) = 1;
Data.Cost.piQMtras(v0) = 1; 

Data.Cost.rhopWi = zeros(length(Data.Red.Branch.T),1);
Data.Cost.rhomqWi = zeros(length(Data.Red.Branch.T),1);
Data.Cost.rhoMqWi = zeros(length(Data.Red.Branch.T),1);
Data.Cost.cdv = ones(length(Data.Red.Branch.T),1);

Data.Cost.piPTras = Data.Cost.piPTras * rhoP_ct * piPTrasHor';
Data.Cost.piQmtras = Data.Cost.piQmtras * rhoQ_ct * piPTrasHor';
Data.Cost.piQMtras = Data.Cost.piQMtras * rhoQ_ct * piPTrasHor';

Data.Cost.m = m*mHor';
Data.Cost.delta = delta;
