function [Data] = cargaEolicosDefault(Data)

	n = size(Data.Red.Branch.T, 1);
    et = size(Data.Red.Bus.pCLow,2);

	Data.Gen.DFIG.Tg = zeros(5, 5);
	Data.Gen.DFIG.Tg(1, 2) = 1;
	Data.Gen.DFIG.Tg(1, 3) = 1;
	Data.Gen.DFIG.Tg(4, 5) = 1;

	Data.Gen.DFIG.I = zeros(n, 1);
    
	Data.Gen.DFIG.C_plb = zeros(n, 1);
	Data.Gen.DFIG.G = zeros(n, 1);
	Data.Gen.DFIG.N_er = zeros(n, 1);
	Data.Gen.DFIG.Np = zeros(n, 1);
	Data.Gen.DFIG.Omega = zeros(n, 1);
	Data.Gen.DFIG.PQnormIE = zeros(n, et);
	Data.Gen.DFIG.PQnormIF = zeros(n, et);
	Data.Gen.DFIG.P_nMec = zeros(n, 1);
	Data.Gen.DFIG.R_ = zeros(n, 1);
	Data.Gen.DFIG.c1 = zeros(n, 1);
	Data.Gen.DFIG.c2 = zeros(n, 1);
	Data.Gen.DFIG.c3 = zeros(n, 1);
	Data.Gen.DFIG.c4 = zeros(n, 1);
	Data.Gen.DFIG.c5 = zeros(n, 1);
	Data.Gen.DFIG.c6 = zeros(n, 1);
	Data.Gen.DFIG.c7 = zeros(n, 1);
	Data.Gen.DFIG.crF = zeros(n, et);
	Data.Gen.DFIG.crR = zeros(n, et);
	Data.Gen.DFIG.cvF = zeros(n, et);
	Data.Gen.DFIG.cvR = zeros(n, et);
	Data.Gen.DFIG.lTopIE = zeros(n, et);
	Data.Gen.DFIG.lTopIF = zeros(n, et);
	Data.Gen.DFIG.lTopOR = zeros(n, et);
	Data.Gen.DFIG.lambda_opt = zeros(n, 1);
	Data.Gen.DFIG.rIE = zeros(n, et);
	Data.Gen.DFIG.rIF = zeros(n, et);
	Data.Gen.DFIG.rOR = zeros(n, et);
	Data.Gen.DFIG.rho = zeros(n, 1);
	Data.Gen.DFIG.sTopF = zeros(n, et);
	Data.Gen.DFIG.sTopR = zeros(n, et);
	Data.Gen.DFIG.uLowE = zeros(n, et);
	Data.Gen.DFIG.uLowF = zeros(n, et);
	Data.Gen.DFIG.uTopE = zeros(n, et);
	Data.Gen.DFIG.uTopF = zeros(n, et);
	Data.Gen.DFIG.vM = zeros(n, 1);
	Data.Gen.DFIG.vMM = zeros(n, 1);
	Data.Gen.DFIG.vm = zeros(n, 1);
	Data.Gen.DFIG.vmm = zeros(n, 1);
	Data.Gen.DFIG.vv = zeros(n, 1);
	Data.Gen.DFIG.w = zeros(n, 1);
	Data.Gen.DFIG.ws = zeros(n, 1);
	Data.Gen.DFIG.xIE = zeros(n, et);
	Data.Gen.DFIG.xIF = zeros(n, et);
	Data.Gen.DFIG.xOR = zeros(n, et);
	Data.Gen.DFIG.xiTopF = zeros(n, et);
	Data.Gen.DFIG.xiTopR = zeros(n, et);

	Data.Gen.DFIG.P_mec = zeros(n, et);
	Data.Gen.DFIG.n_ = zeros(n, et);
    