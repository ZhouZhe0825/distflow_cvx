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
	DataM.Red.Branch.NtrLow = repmat(full(DataM.Red.Branch.NtrLow), [1 1 lenEt]);
	DataM.Red.Branch.NtrTop = repmat(full(DataM.Red.Branch.NtrTop), [1 1 lenEt]);
	DataM.Red.Branch.Tap = repmat(full(DataM.Red.Branch.Tap), [1 1 lenEt]);

	DataM.Red.Branch.lTop = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.lTop);
	DataM.Red.Branch.r = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.r);
	DataM.Red.Branch.x = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.x);
	DataM.Red.Branch.yLow = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.yLow);
	DataM.Red.Branch.yTop = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.yTop);
	DataM.Red.Branch.NtrLow = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.NtrLow);
	DataM.Red.Branch.NtrTop = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.NtrTop);
	DataM.Red.Branch.NtrIni = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.NtrIni);
	DataM.Red.Branch.Tap = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.Tap);
    
    
	a1 = repmat(DataM.Red.Bus.alpha(:,1), [1 lenEt 2]);
	a1(:,:,2) = repmat(DataM.Red.Bus.alpha(:,2), [1 lenEt]);
	DataM.Red.Bus.alpha = a1;
	DataM.Red.Bus.NcpLow = DataM.Red.Bus.NcpLow * one;
	DataM.Red.Bus.NcpTop = DataM.Red.Bus.NcpTop * one;
	DataM.Red.Bus.Cap = DataM.Red.Bus.Cap * one;
	DataM.Red.Bus.pCLow = full(DataM.Red.Bus.pCLow(:,Et));
	DataM.Red.Bus.qCLow = full(DataM.Red.Bus.qCLow(:,Et));
	DataM.Red.Bus.uLow = DataM.Red.Bus.uLow * one;
	DataM.Red.Bus.uTop = DataM.Red.Bus.uTop * one;
	
	DataM.ClNI.pC = DataM.ClNI.pC * one;
	DataM.ClNI.qC = DataM.ClNI.qC * one;

	
	
	
	DataM.Cost.cdv = DataM.Cost.cdv(:,Et);
	DataM.Cost.cCap = DataM.Cost.cCap * one;
	DataM.Cost.cBas = DataM.Cost.cBas(:,Et);
	DataM.Cost.cTap = repmat(full(DataM.Cost.cTap), [1 1 lenEt]);
	DataM.Cost.cTap = NxNxT2MxT(VertI,VertJ,DataM.Cost.cTap);
	DataM.Cost.cY = repmat(full(DataM.Cost.cY), [1 1 lenEt]);
	DataM.Cost.cY = NxNxT2MxT(VertI,VertJ,DataM.Cost.cY);
    DataM.Cost.delta = DataM.Cost.delta(:,Et);
	DataM.Cost.m = DataM.Cost.m(:,Et);
	DataM.Cost.piPTras = DataM.Cost.piPTras(:,Et);
	DataM.Cost.piQmtras = DataM.Cost.piQmtras(:,Et);
	DataM.Cost.piQMtras = DataM.Cost.piQMtras(:,Et);
	DataM.Cost.rhopPv	 = DataM.Cost.rhopPv(:,Et);
	DataM.Cost.rhomqPv	 = DataM.Cost.rhomqPv(:,Et);
	DataM.Cost.rhoMqPv	 = DataM.Cost.rhoMqPv(:,Et);
	DataM.Cost.rhopWi = DataM.Cost.rhopWi(:,Et);
	DataM.Cost.rhomqWi = DataM.Cost.rhomqWi(:,Et);
	DataM.Cost.rhoMqWi = DataM.Cost.rhoMqWi(:,Et);




	if lenWn > 0
		DataM.Gen.DFIG.crF = DataM.Gen.DFIG.crF(:,Et);
		DataM.Gen.DFIG.crR = DataM.Gen.DFIG.crR(:,Et);
		DataM.Gen.DFIG.cvF = DataM.Gen.DFIG.cvF(:,Et);
		DataM.Gen.DFIG.cvR = DataM.Gen.DFIG.cvR(:,Et);
		DataM.Gen.DFIG.lTopIE = DataM.Gen.DFIG.lTopIE(:,Et);
		DataM.Gen.DFIG.lTopIF = DataM.Gen.DFIG.lTopIF(:,Et);
		DataM.Gen.DFIG.lTopOR = DataM.Gen.DFIG.lTopOR(:,Et);
		DataM.Gen.DFIG.n_	 = 	DataM.Gen.DFIG.n_(:,Et);
		DataM.Gen.DFIG.P_mec	 = 	DataM.Gen.DFIG.P_mec(:,Et);
		DataM.Gen.DFIG.PQnormIE = DataM.Gen.DFIG.PQnormIE(:,Et);
		DataM.Gen.DFIG.PQnormIF = DataM.Gen.DFIG.PQnormIF(:,Et);
		DataM.Gen.DFIG.rIE = DataM.Gen.DFIG.rIE(:,Et);
		DataM.Gen.DFIG.rIF = DataM.Gen.DFIG.rIF(:,Et);
		DataM.Gen.DFIG.rOR = DataM.Gen.DFIG.rOR(:,Et);
		DataM.Gen.DFIG.sTopF = DataM.Gen.DFIG.sTopF(:,Et);
		DataM.Gen.DFIG.sTopR = DataM.Gen.DFIG.sTopR(:,Et);
		DataM.Gen.DFIG.xIE = DataM.Gen.DFIG.xIE(:,Et);
		DataM.Gen.DFIG.xIF = DataM.Gen.DFIG.xIF(:,Et);
		DataM.Gen.DFIG.xOR = DataM.Gen.DFIG.xOR(:,Et);
		DataM.Gen.DFIG.xiTopF = DataM.Gen.DFIG.xiTopF(:,Et);
		DataM.Gen.DFIG.xiTopR = DataM.Gen.DFIG.xiTopR(:,Et);
		DataM.Gen.DFIG.uLowE = DataM.Gen.DFIG.uLowE(:,Et);
		DataM.Gen.DFIG.uLowF = DataM.Gen.DFIG.uLowF(:,Et);
		DataM.Gen.DFIG.uTopE = DataM.Gen.DFIG.uTopE(:,Et);
		DataM.Gen.DFIG.uTopF = DataM.Gen.DFIG.uTopF(:,Et);
		DataM.Gen.DFIG.qWiLow = DataM.Gen.DFIG.qWiLow(:,Et);
		DataM.Gen.DFIG.qWiTop = DataM.Gen.DFIG.qWiTop(:,Et);
	end

	DataM.Gen.Pv.cr = DataM.Gen.Pv.cr * one;
	DataM.Gen.Pv.cv = DataM.Gen.Pv.cv * one;
	DataM.Gen.Pv.pPvg = repmat(DataM.Gen.Pv.pPvg(:,Et), [size(VertI,2) 1]);
	DataM.Gen.Pv.sTop = DataM.Gen.Pv.sTop * one;
	DataM.Gen.Pv.xiTop = DataM.Gen.Pv.xiTop * one;

	DataM.Gen.Tras.pgLow = DataM.Gen.Tras.pgLow(:,Et);
	DataM.Gen.Tras.pgTop = DataM.Gen.Tras.pgTop(:,Et);
	DataM.Gen.Tras.qgLow = DataM.Gen.Tras.qgLow(:,Et);
	DataM.Gen.Tras.qgTop = DataM.Gen.Tras.qgTop(:,Et);

	DataM.Gen.Basic.pgLow = DataM.Gen.Basic.pgLow(:,Et);
	DataM.Gen.Basic.pgTop = DataM.Gen.Basic.pgTop(:,Et);
	DataM.Gen.Basic.qgLow = DataM.Gen.Basic.qgLow(:,Et);
	DataM.Gen.Basic.qgTop = DataM.Gen.Basic.qgTop(:,Et);
	
	DataM.St.AC.beta = DataM.St.AC.beta(:,Et);
	DataM.St.AC.epsilon = DataM.St.AC.epsilon(:,Et);
	DataM.St.AC.eta = DataM.St.AC.eta(:,Et);
	DataM.St.AC.tempLow = DataM.St.AC.tempLow(:,Et);
	DataM.St.AC.tempPref = DataM.St.AC.tempPref(:,Et);
	DataM.St.AC.tempTop = DataM.St.AC.tempTop(:,Et);

	DataM.St.Bat.beta = DataM.St.Bat.beta * one;
	DataM.St.Bat.cr = DataM.St.Bat.cr * one;
	DataM.St.Bat.cv = DataM.St.Bat.cv * one;
	DataM.St.Bat.ELow = DataM.St.Bat.ELow * one;
	DataM.St.Bat.ETop = DataM.St.Bat.ETop * one;
	DataM.St.Bat.EPref = DataM.St.Bat.EPref * one;
	DataM.St.Bat.epsilon = DataM.St.Bat.epsilon * one;
	DataM.St.Bat.etaC = DataM.St.Bat.etaC * one;
	DataM.St.Bat.etaD = DataM.St.Bat.etaD * one;
	DataM.St.Bat.gama = DataM.St.Bat.gama * one;
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

	DataM.temp = repmat(DataM.temp(:,Et), [size(VertI,2) 1]);

	DataM.Util.betaE = DataM.Util.betaE(:,Et);
	DataM.Util.betaT = DataM.Util.betaT(:,Et);
	DataM.Util.pzCnLow = DataM.Util.pzCnLow(:,Et,:);
	DataM.Util.pzCnLowE = DataM.Util.pzCnLowE * one;
	DataM.Util.pzCnPref = DataM.Util.pzCnPref(:,Et,:);
	DataM.Util.pzCnPrefE = DataM.Util.pzCnPrefE * one;
	DataM.Util.pzCnTop = DataM.Util.pzCnTop(:,Et,:);
	DataM.Util.pzCnTopE = DataM.Util.pzCnTopE * one;
	DataM.Util.qzCnLowE = DataM.Util.qzCnLowE * one;
	DataM.Util.qzCnTopE = DataM.Util.qzCnTopE * one;

end