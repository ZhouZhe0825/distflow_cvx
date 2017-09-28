function [Data] = DataDefault(n,et,app)

	Data.ClNI.I = zeros(n,1);
	Data.ClNI.d = zeros(n,1);
	Data.ClNI.nMultipLow = zeros(n,1);
	Data.ClNI.nMultipTop = zeros(n,1);
	Data.ClNI.pC = zeros(n,1);
	Data.ClNI.qC = zeros(n,1);
	Data.Cost.cCap = zeros(n,1);
	Data.Cost.cTap = zeros(n,n);
	Data.Cost.cBas = zeros(n,et);
	Data.Cost.cY = zeros(n,n);
	Data.Cost.cdv = zeros(n,et);
	Data.Cost.delta = zeros(n,et);
	Data.Cost.m = zeros(n,et);
	Data.Cost.piPTras = zeros(n,et);
	Data.Cost.piQMtras = zeros(n,et);
	Data.Cost.piQmtras = zeros(n,et);
	Data.Cost.rhoMqPv = zeros(n,et);
	Data.Cost.rhoMqWi = zeros(n,et);
	Data.Cost.rhomqPv = zeros(n,et);
	Data.Cost.rhomqWi = zeros(n,et);
	Data.Cost.rhopPv = zeros(n,et);
	Data.Cost.rhopWi = zeros(n,et);
    Data.Gen.Basic.pgLow = zeros(n,et);
    Data.Gen.Basic.qgLow = zeros(n,et);
    Data.Gen.Basic.pgTop = zeros(n,et);
    Data.Gen.Basic.qgTop = zeros(n,et);
    Data.Gen.Basic.pgIni = zeros(n,1);
    Data.Gen.Basic.qgIni = zeros(n,1);
    Data.Gen.Basic.I = zeros(n,1);
    Data.Gen.DFIG.C_plb = zeros(n,1);
	Data.Gen.DFIG.G = zeros(n,1);
	Data.Gen.DFIG.I = zeros(n,1);
	Data.Gen.DFIG.N_er = zeros(n,1);
	Data.Gen.DFIG.Np = zeros(n,1);
	Data.Gen.DFIG.Omega = zeros(n,1);
	Data.Gen.DFIG.PQnormIE = zeros(n,et);
	Data.Gen.DFIG.PQnormIF = zeros(n,et);
	Data.Gen.DFIG.P_mec = zeros(n,et);
	Data.Gen.DFIG.P_nMec = zeros(n,1);
	Data.Gen.DFIG.R_ = zeros(n,1);
	Data.Gen.DFIG.c1 = zeros(n,1);
	Data.Gen.DFIG.c2 = zeros(n,1);
	Data.Gen.DFIG.c3 = zeros(n,1);
	Data.Gen.DFIG.c4 = zeros(n,1);
	Data.Gen.DFIG.c5 = zeros(n,1);
	Data.Gen.DFIG.c6 = zeros(n,1);
	Data.Gen.DFIG.c7 = zeros(n,1);
	Data.Gen.DFIG.crF = zeros(n,et);
	Data.Gen.DFIG.crR = zeros(n,et);
	Data.Gen.DFIG.cvF = zeros(n,et);
	Data.Gen.DFIG.cvR = zeros(n,et);
	Data.Gen.DFIG.lTopIE = zeros(n,et);
	Data.Gen.DFIG.lTopIF = zeros(n,et);
	Data.Gen.DFIG.lTopOR = zeros(n,et);
	Data.Gen.DFIG.lambda_opt = zeros(n,1);
	Data.Gen.DFIG.n_ = zeros(n,et);
	Data.Gen.DFIG.rIE = zeros(n,et);
	Data.Gen.DFIG.rIF = zeros(n,et);
	Data.Gen.DFIG.rOR = zeros(n,et);
	Data.Gen.DFIG.rho = zeros(n,1);
	Data.Gen.DFIG.sTopF = zeros(n,et);
	Data.Gen.DFIG.sTopR = zeros(n,et);
	Data.Gen.DFIG.uLowE = zeros(n,et);
	Data.Gen.DFIG.uLowF = zeros(n,et);
	Data.Gen.DFIG.uTopE = zeros(n,et);
	Data.Gen.DFIG.uTopF = zeros(n,et);
	Data.Gen.DFIG.vM = zeros(n,1);
	Data.Gen.DFIG.vMM = zeros(n,1);
	Data.Gen.DFIG.vm = zeros(n,1);
	Data.Gen.DFIG.vmm = zeros(n,1);
	Data.Gen.DFIG.vv = zeros(n,1);
	Data.Gen.DFIG.w = zeros(n,1);
	Data.Gen.DFIG.ws = zeros(n,1);
	Data.Gen.DFIG.xIE = zeros(n,et);
	Data.Gen.DFIG.xIF = zeros(n,et);
	Data.Gen.DFIG.xOR = zeros(n,et);
	Data.Gen.DFIG.xiTopF = zeros(n,et);
	Data.Gen.DFIG.xiTopR = zeros(n,et);
    Data.Gen.DFIG.pgIni = zeros(n,1);
    Data.Gen.DFIG.qgIni = zeros(n,1);
	Data.Gen.DFIG.qWiLow = zeros(n,et);
	Data.Gen.DFIG.qWiTop = zeros(n,et);
	Data.Gen.Pv.I = zeros(n,1);
	Data.Gen.Pv.cr = zeros(n,1);
	Data.Gen.Pv.cv = zeros(n,1);
	Data.Gen.Pv.pPvg = zeros(n,et);
	Data.Gen.Pv.pgTop = zeros(n,1);
	Data.Gen.Pv.sTop = zeros(n,1);
	Data.Gen.Pv.xiTop = zeros(n,1);
    Data.Gen.Pv.pgIni = zeros(n,1);
    Data.Gen.Pv.qgIni = zeros(n,1);
    Data.Gen.Tras.I = zeros(n,1);
	Data.Gen.Tras.pgLow = zeros(n,et);
	Data.Gen.Tras.pgTop = zeros(n,et);
	Data.Gen.Tras.qgLow = zeros(n,et);
	Data.Gen.Tras.qgTop = zeros(n,et);
    Data.Gen.Tras.pgIni = zeros(n,1);
    Data.Gen.Tras.qgIni = zeros(n,1);
	Data.Red.Branch.T = zeros(n,n);
	Data.Red.Branch.Tswitches = zeros(n,n);
	Data.Red.Branch.lTop = zeros(n,n);
	Data.Red.Branch.r = zeros(n,n);
	Data.Red.Branch.x = zeros(n,n);
	Data.Red.Branch.yLow = zeros(n,n);
	Data.Red.Branch.yTop = zeros(n,n);
	Data.Red.Branch.Itap = zeros(n,n);
	Data.Red.Branch.Itreg = zeros(n,n);
	Data.Red.Branch.NtrIni = zeros(n,n);
	Data.Red.Branch.NtrLow = zeros(n,n);
	Data.Red.Branch.NtrTop = zeros(n,n);
	Data.Red.Branch.Tap = zeros(n,n);
	Data.Red.Bus.Cap = zeros(n,1);
	Data.Red.Bus.Icap = zeros(n,1);
	Data.Red.Bus.Icons = zeros(n,1);
    Data.Red.Bus.NcpIni = zeros(n,1);
	Data.Red.Bus.NcpLow = zeros(n,1);
	Data.Red.Bus.NcpTop = zeros(n,1);
	Data.Red.Bus.alpha = zeros(n,app);
	Data.Red.Bus.pCLow = zeros(n,et);
	Data.Red.Bus.qCLow = zeros(n,et);
	Data.Red.Bus.uLow = zeros(n,1);
	Data.Red.Bus.uTop = zeros(n,1);
	Data.St.AC.I = zeros(n,1);
	Data.St.AC.beta = zeros(n,et);
	Data.St.AC.epsilon = zeros(n,et);
	Data.St.AC.eta = zeros(n,et);
	Data.St.AC.tempIni = zeros(n,1);
	Data.St.AC.tempLow = zeros(n,et);
	Data.St.AC.tempPref = zeros(n,et);
	Data.St.AC.tempTop = zeros(n,et);
	Data.St.Bat.EIni = zeros(n,1);
	Data.St.Bat.ELow = zeros(n,1);
	Data.St.Bat.ETop = zeros(n,1);
	Data.St.Bat.EPref = zeros(n,1);
	Data.St.Bat.I = zeros(n,1);
	Data.St.Bat.beta = zeros(n,1);
	Data.St.Bat.cr = zeros(n,1);
	Data.St.Bat.cv = zeros(n,1);
	Data.St.Bat.epsilon = zeros(n,1);
	Data.St.Bat.etaC = zeros(n,1);
	Data.St.Bat.etaD = zeros(n,1);
	Data.St.Bat.gama = zeros(n,1);
	Data.St.Bat.kapa = zeros(n,1);
	Data.St.Bat.m1 = zeros(n,1);
	Data.St.Bat.m2 = zeros(n,1);
	Data.St.Bat.m3 = zeros(n,1);
	Data.St.Bat.pgLow = zeros(n,1);
	Data.St.Bat.pgTop = zeros(n,1);
	Data.St.Bat.sTop = zeros(n,1);
	Data.St.Bat.wOm = zeros(n,1);
	Data.St.Bat.wU = zeros(n,1);
	Data.St.Bat.xiTop = zeros(n,1);
    Data.St.Bat.pgIni = zeros(n,1);
    Data.St.Bat.qgIni = zeros(n,1);
	Data.Util.betaE = zeros(n,et);
	Data.Util.betaT = zeros(n,et);
	Data.Util.pzCnLow = zeros(n,et,app);
	Data.Util.pzCnLowE = zeros(n,1);
	Data.Util.pzCnPref = zeros(n,et,app);
	Data.Util.pzCnPrefE = zeros(n,1);
	Data.Util.pzCnTop = zeros(n,et,app);
	Data.Util.pzCnTopE = zeros(n,1);
	Data.Util.qzCnLowE = zeros(n,1);
	Data.Util.qzCnPrefE = zeros(n,1);
	Data.Util.qzCnTopE = zeros(n,1);
	Data.Util.tgPhi = zeros(app,1);
	Data.dt = .25;
	Data.temp = zeros(1,et);

end