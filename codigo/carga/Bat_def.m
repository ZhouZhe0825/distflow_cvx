function [Data] = Bat_def(Data,EIni,batNodes)

	cv_ct = 0.001;
	cr_ct = 0.001;
	epsilon_ct = 0.1;
	eta_ct = 1;
	pgTop_ct = 0.5;
	pgLow_ct = -0.5;
	sTop_ct = 0.095;
	ETop_ct = 2;
	ELow_ct = 0;
	beta_ct = 0.2;
	wU_ct = 1;
	wOm_ct = 1;
	m1_ct = 1;
	m2_ct = 0.75;
	m3_ct = 0.5;
	kapa_ct = 0.2;
	gama_ct = 0.6;


	Data.St.Bat.I = zeros(size(Data.Red.Branch.T,1),1);
	Data.St.Bat.I(batNodes) = 1;
	Data.St.Bat.cv = Data.St.Bat.I * cv_ct; % temperatura minima, por nodo
	Data.St.Bat.cr = Data.St.Bat.I * cr_ct; % temperatura maxima, por nodo
	Data.St.Bat.epsilon = Data.St.Bat.I * epsilon_ct; % epsilon por nodo
	Data.St.Bat.eta = Data.St.Bat.I * eta_ct; % eta por nodo
	Data.St.Bat.pgTop = Data.St.Bat.I * pgTop_ct; % temperatura minima, por nodo
	Data.St.Bat.pgLow = Data.St.Bat.I * pgLow_ct; % temperatura maxima, por nodo
	Data.St.Bat.sTop = Data.St.Bat.I * sTop_ct; % temperatura minima, por nodo
	Data.St.Bat.ETop = Data.St.Bat.I * ETop_ct; % temperatura minima, por nodo
	Data.St.Bat.ELow = Data.St.Bat.I * ELow_ct; % temperatura maxima, por nodo
	Data.St.Bat.beta = Data.St.Bat.I * beta_ct; %betaAC por nodo
	Data.St.Bat.wU = Data.St.Bat.I * wU_ct; %betaAC por nodo
	Data.St.Bat.wOm = Data.St.Bat.I * wOm_ct; %aAC por nodo
	Data.St.Bat.m1 = Data.St.Bat.I * m1_ct; %aAC por nodo
	Data.St.Bat.m2 = Data.St.Bat.I * m2_ct; %aAC por nodo
	Data.St.Bat.m3 = Data.St.Bat.I * m3_ct; %aAC por nodo

	Data.St.Bat.xiTop = Data.St.Bat.sTop.^2; % temperatura minima, por nodo
	
	Data.St.Bat.kapa = kapa_ct; %aAC por nodo
	Data.St.Bat.gama = gama_ct; %aAC por nodo

	Data.St.Bat.EIni = Data.St.Bat.I;
	Data.St.Bat.EIni(batNodes) = EIni;

	
end