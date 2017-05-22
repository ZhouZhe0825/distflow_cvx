function [DataM] = reshapeDataM(Data, Config)

	VertI = VertIMat(Data.Red.Branch.T);
	VertJ = VertJMat(Data.Red.Branch.T);
	iniEt = Config.iniEtapa;
	finEt = iniEt+Config.Etapas-1;
	Et = (iniEt:finEt);
	lenEt = length(Et);
	one = ones(1,lenEt);

	lenWn = length(find(Data.Gen.DFIG.I == 1));

	DataM = Data;

	
	
	DataM.Gen.Tras.pgLow = full(DataM.Gen.Tras.pgLow);
	DataM.Gen.Tras.qgLow = full(DataM.Gen.Tras.qgLow);
	DataM.Gen.Tras.pgTop = full(DataM.Gen.Tras.pgTop);
	DataM.Gen.Tras.qgTop = full(DataM.Gen.Tras.qgTop);


	DataM.Red.Branch.lTop = repmat(full(DataM.Red.Branch.lTop), [1 1 lenEt]);
	DataM.Red.Branch.r = repmat(full(DataM.Red.Branch.r), [1 1 lenEt]);
	DataM.Red.Branch.x = repmat(full(DataM.Red.Branch.x), [1 1 lenEt]);
	DataM.Red.Branch.yLow = repmat(full(DataM.Red.Branch.yLow), [1 1 lenEt]);
	DataM.Red.Branch.yTop = repmat(full(DataM.Red.Branch.yTop), [1 1 lenEt]);

	DataM.Red.Branch.lTop = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.lTop);
	DataM.Red.Branch.r = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.r);
	DataM.Red.Branch.x = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.x);
	DataM.Red.Branch.yLow = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.yLow);
	DataM.Red.Branch.yTop = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.yTop);

	a1 = repmat(DataM.Red.Bus.alpha(:,1), [1 lenEt 2]);
	a1(:,:,2) = repmat(DataM.Red.Bus.alpha(:,2), [1 lenEt]);
	DataM.Red.Bus.alpha = a1;
	DataM.Red.Bus.CapLow = DataM.Red.Bus.CapLow * one;
	DataM.Red.Bus.CapTop = DataM.Red.Bus.CapTop * one;
	DataM.Red.Bus.indCap = DataM.Red.Bus.indCap * one;
	DataM.Red.Bus.indTap = DataM.Red.Bus.indTap * one;
	DataM.Red.Bus.Ncp = DataM.Red.Bus.Ncp * one;
	DataM.Red.Bus.Ntr = DataM.Red.Bus.Ntr * one;
	DataM.Red.Bus.pCLow = full(DataM.Red.Bus.pCLow(:,(1:lenEt)));
	DataM.Red.Bus.qCLow = full(DataM.Red.Bus.qCLow(:,(1:lenEt)));
	DataM.Red.Bus.TapLow = DataM.Red.Bus.TapLow * one;
	DataM.Red.Bus.TapTop = DataM.Red.Bus.TapTop * one;
	DataM.Red.Bus.uLow = DataM.Red.Bus.uLow * one;
	DataM.Red.Bus.uTop = DataM.Red.Bus.uTop * one;

	
	DataM.ClNI.pC = DataM.ClNI.pC * one;
	DataM.ClNI.qC = DataM.ClNI.qC * one;

	
	
	
	DataM.Cost.cdv = DataM.Cost.cdv(:,(1:lenEt));
	DataM.Cost.cCap = DataM.Cost.cCap * one;
	DataM.Cost.cTap = DataM.Cost.cTap * one;
	DataM.Cost.cY = repmat(full(DataM.Cost.cY), [1 1 lenEt]);
	DataM.Cost.cY = NxNxT2MxT(VertI,VertJ,DataM.Cost.cY);
    DataM.Cost.delta = DataM.Cost.delta(:,(1:lenEt));
	DataM.Cost.m = DataM.Cost.m(:,(1:lenEt));
	DataM.Cost.piPTras = DataM.Cost.piPTras(:,(1:lenEt));
	DataM.Cost.piQmtras = DataM.Cost.piQmtras(:,(1:lenEt));
	DataM.Cost.piQMtras = DataM.Cost.piQMtras(:,(1:lenEt));
	DataM.Cost.rhopPv	 = DataM.Cost.rhopPv(:,(1:lenEt));
	DataM.Cost.rhomqPv	 = DataM.Cost.rhomqPv(:,(1:lenEt));
	DataM.Cost.rhoMqPv	 = DataM.Cost.rhoMqPv(:,(1:lenEt));
	DataM.Cost.rhopWi = DataM.Cost.rhopWi(:,(1:lenEt));
	DataM.Cost.rhomqWi = DataM.Cost.rhomqWi(:,(1:lenEt));
	DataM.Cost.rhoMqWi = DataM.Cost.rhoMqWi(:,(1:lenEt));




	if lenWn > 0
		DataM.Gen.DFIG.crF = DataM.Gen.DFIG.crF(:,(1:lenEt));
		DataM.Gen.DFIG.crR = DataM.Gen.DFIG.crR(:,(1:lenEt));
		DataM.Gen.DFIG.cvF = DataM.Gen.DFIG.cvF(:,(1:lenEt));
		DataM.Gen.DFIG.cvR = DataM.Gen.DFIG.cvR(:,(1:lenEt));
		DataM.Gen.DFIG.lTopIE = DataM.Gen.DFIG.lTopIE(:,(1:lenEt));
		DataM.Gen.DFIG.lTopIF = DataM.Gen.DFIG.lTopIF(:,(1:lenEt));
		DataM.Gen.DFIG.lTopOR = DataM.Gen.DFIG.lTopOR(:,(1:lenEt));
		DataM.Gen.DFIG.n_	 = 	DataM.Gen.DFIG.n_(:,(1:lenEt));
		DataM.Gen.DFIG.P_mec	 = 	DataM.Gen.DFIG.P_mec(:,(1:lenEt));
		DataM.Gen.DFIG.PQnormIE = DataM.Gen.DFIG.PQnormIE(:,(1:lenEt));
		DataM.Gen.DFIG.PQnormIF = DataM.Gen.DFIG.PQnormIF(:,(1:lenEt));
		DataM.Gen.DFIG.rIE = DataM.Gen.DFIG.rIE(:,(1:lenEt));
		DataM.Gen.DFIG.rIF = DataM.Gen.DFIG.rIF(:,(1:lenEt));
		DataM.Gen.DFIG.rOR = DataM.Gen.DFIG.rOR(:,(1:lenEt));
		DataM.Gen.DFIG.sTopF = DataM.Gen.DFIG.sTopF(:,(1:lenEt));
		DataM.Gen.DFIG.sTopR = DataM.Gen.DFIG.sTopR(:,(1:lenEt));
		DataM.Gen.DFIG.xIE = DataM.Gen.DFIG.xIE(:,(1:lenEt));
		DataM.Gen.DFIG.xIF = DataM.Gen.DFIG.xIF(:,(1:lenEt));
		DataM.Gen.DFIG.xOR = DataM.Gen.DFIG.xOR(:,(1:lenEt));
		DataM.Gen.DFIG.xiTopF = DataM.Gen.DFIG.xiTopF(:,(1:lenEt));
		DataM.Gen.DFIG.xiTopR = DataM.Gen.DFIG.xiTopR(:,(1:lenEt));
		DataM.Gen.DFIG.uLowE = DataM.Gen.DFIG.uLowE(:,(1:lenEt));
		DataM.Gen.DFIG.uLowF = DataM.Gen.DFIG.uLowF(:,(1:lenEt));
		DataM.Gen.DFIG.uTopE = DataM.Gen.DFIG.uTopE(:,(1:lenEt));
		DataM.Gen.DFIG.uTopF = DataM.Gen.DFIG.uTopF(:,(1:lenEt));
	end

	DataM.Gen.Pv.cr = DataM.Gen.Pv.cr * one;
	DataM.Gen.Pv.cv = DataM.Gen.Pv.cv * one;
	DataM.Gen.Pv.I = DataM.Gen.Pv.I * one;
	DataM.Gen.Pv.pPvg = repmat(DataM.Gen.Pv.pPvg(:,(1:lenEt)), [size(VertI,2) 1]);
	DataM.Gen.Pv.sTop = DataM.Gen.Pv.sTop * one;
	DataM.Gen.Pv.xiTop = DataM.Gen.Pv.xiTop * one;

	DataM.Gen.Tras.pgLow = DataM.Gen.Tras.pgLow * one;
	DataM.Gen.Tras.pgTop = DataM.Gen.Tras.pgTop * one;
	DataM.Gen.Tras.qgLow = DataM.Gen.Tras.qgLow * one;
	DataM.Gen.Tras.qgTop = DataM.Gen.Tras.qgTop * one;

	
	
	
	DataM.St.AC.beta = DataM.St.AC.beta(:,(1:lenEt));
	DataM.St.AC.epsilon = DataM.St.AC.epsilon(:,(1:lenEt));
	DataM.St.AC.eta = DataM.St.AC.eta(:,(1:lenEt));
	DataM.St.AC.tempLow = DataM.St.AC.tempLow(:,(1:lenEt));
	DataM.St.AC.tempPref = DataM.St.AC.tempPref(:,(1:lenEt));
	DataM.St.AC.tempTop = DataM.St.AC.tempTop(:,(1:lenEt));

	DataM.St.Bat.beta = DataM.St.Bat.beta * one;
	DataM.St.Bat.cr = DataM.St.Bat.cr * one;
	DataM.St.Bat.cv = DataM.St.Bat.cv * one;
	DataM.St.Bat.ELow = DataM.St.Bat.ELow * one;
	DataM.St.Bat.ETop = DataM.St.Bat.ETop * one;
	DataM.St.Bat.epsilon = DataM.St.Bat.epsilon * one;
	DataM.St.Bat.eta = DataM.St.Bat.eta * one;
	DataM.St.Bat.gama = DataM.St.Bat.gama * one;
	DataM.St.Bat.I = DataM.St.Bat.I * one;
	DataM.St.Bat.kapa = DataM.St.Bat.kapa * one;
	DataM.St.Bat.m1 = DataM.St.Bat.m1 * one;
	DataM.St.Bat.m2 = DataM.St.Bat.m2 * one;
	DataM.St.Bat.m3 = DataM.St.Bat.m3 * one;
	DataM.St.Bat.pgLow = DataM.St.Bat.pgLow * one;
	DataM.St.Bat.pgTop = DataM.St.Bat.pgTop * one;
	DataM.St.Bat.sTop = DataM.St.Bat.sTop * one;
	DataM.St.Bat.xiTop = DataM.St.Bat.xiTop * one;
	DataM.St.Bat.wOm = DataM.St.Bat.wOm * one;
	DataM.St.Bat.wU = DataM.St.Bat.wU * one;

	DataM.temp = repmat(DataM.temp(:,(1:lenEt)), [size(VertI,2) 1]);

	DataM.Util.betaE = DataM.Util.betaE(:,(1:lenEt));
	DataM.Util.betaT = DataM.Util.betaT(:,(1:lenEt));
	DataM.Util.pzCnLow = DataM.Util.pzCnLow(:,(1:lenEt),:);
	DataM.Util.pzCnLowE = DataM.Util.pzCnLowE * one;
	DataM.Util.pzCnPref = DataM.Util.pzCnPref(:,(1:lenEt),:);
	DataM.Util.pzCnPrefE = DataM.Util.pzCnPrefE * one;
	DataM.Util.pzCnTop = DataM.Util.pzCnTop(:,(1:lenEt),:);
	DataM.Util.pzCnTopE = DataM.Util.pzCnTopE * one;
	DataM.Util.qzCnLowE = DataM.Util.qzCnLowE * one;
	DataM.Util.qzCnTopE = DataM.Util.qzCnTopE * one;

end