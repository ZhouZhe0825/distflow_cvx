function [Data] = cargaPvDefault(Data)

	Data.Gen.Pv.I = zeros(size(Data.Gen.Tras.I));
	Data.Gen.Pv.pPvg = zeros(size(Data.Gen.Pv.I,1),size(Data.Red.Bus.pCLow,2));
	Data.Gen.Pv.sTop = Data.Gen.Pv.I;
	Data.Gen.Pv.xiTop = Data.Gen.Pv.I;											
	Data.Gen.Pv.pgTop = Data.Gen.Pv.I;
	Data.Gen.Pv.cv = Data.Gen.Pv.I;
	Data.Gen.Pv.cr = Data.Gen.Pv.I;

end