function [Data] = PvGen_sm(Data,pvNodes)

	sTop_ct = 0.8;
	pgTop_ct = 0.0075;
	cv_ct = .1;
	cr_ct = .1;
	
	pPvg_top = 0.3;
	pPvg_low = .05;

	Data.Gen.Pv.I = zeros(size(Data.Gen.Tras.I));
	Data.Gen.Pv.I(pvNodes) = 1;
	Data.Gen.Pv.pPvg = (Data.temp - min(Data.temp))/max(Data.temp)*pPvg_top + pPvg_low; % potencia de generacion del solar, por sale desde la temperatura
	Data.Gen.Pv.sTop = Data.Gen.Pv.I * sTop_ct;
	Data.Gen.Pv.xiTop = Data.Gen.Pv.sTop .^ 2;											
	Data.Gen.Pv.pgTop = Data.Gen.Pv.I * pgTop_ct;
	Data.Gen.Pv.cv = Data.Gen.Pv.I * cv_ct;
	Data.Gen.Pv.cr = Data.Gen.Pv.I * cr_ct;

end