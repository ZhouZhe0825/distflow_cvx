function [DataM] = reshapeDataM(Data, Config);

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


	DataM.Red.Branch.cY = repmat(full(DataM.Red.Branch.cY), [1 1 lenEt]);
	DataM.Red.Branch.lTop = repmat(full(DataM.Red.Branch.lTop), [1 1 lenEt]);
	DataM.Red.Branch.r = repmat(full(DataM.Red.Branch.r), [1 1 lenEt]);
	DataM.Red.Branch.x = repmat(full(DataM.Red.Branch.x), [1 1 lenEt]);
	DataM.Red.Branch.yLow = repmat(full(DataM.Red.Branch.yLow), [1 1 lenEt]);
	DataM.Red.Branch.yTop = repmat(full(DataM.Red.Branch.yTop), [1 1 lenEt]);

	DataM.Red.Branch.cY = NxNxT2MxT(VertI,VertJ,DataM.Red.Branch.cY);
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

	
	
	
	DataM.Cost.cdv = full(DataM.Cost.cdv) * one;
	DataM.Cost.m = DataM.Cost.m(:,(1:lenEt));
	DataM.Cost.piPTras = DataM.Cost.piPTras(:,(1:lenEt));
	DataM.Cost.piQmtras = DataM.Cost.piQmtras(:,(1:lenEt));
	DataM.Cost.piQMtras = DataM.Cost.piQMtras(:,(1:lenEt));
	DataM.Cost.rhopPv	 = DataM.Cost.rhopPv * one;
	DataM.Cost.rhomqPv	 = DataM.Cost.rhomqPv * one;
	DataM.Cost.rhoMqPv	 = DataM.Cost.rhoMqPv * one;
	DataM.Cost.rhopWi = DataM.Cost.rhopWi * one;
	DataM.Cost.rhomqWi = DataM.Cost.rhomqWi * one;
	DataM.Cost.rhoMqWi = DataM.Cost.rhoMqWi * one;




	if lenWn > 0
		DataM.Gen.DFIG.crF = DataM.Gen.DFIG.cr(3,:)' * one;
		DataM.Gen.DFIG.crR = DataM.Gen.DFIG.cr(5,:)' * one;
		DataM.Gen.DFIG.cvF = DataM.Gen.DFIG.cv(3,:)' * one;
		DataM.Gen.DFIG.cvR = DataM.Gen.DFIG.cv(5,:)' * one;
		DataM.Gen.DFIG.lTopIE = squeeze(DataM.Gen.DFIG.lTop(1,2,:)) * one;
		DataM.Gen.DFIG.lTopIF = squeeze(DataM.Gen.DFIG.lTop(1,3,:)) * one;
		DataM.Gen.DFIG.lTopOR = squeeze(DataM.Gen.DFIG.lTop(4,5,:)) * one;
		DataM.Gen.DFIG.n_	 = 	repmat(DataM.Gen.DFIG.n_(:,Et), [lenWn 1]);
		DataM.Gen.DFIG.P_mec	 = 	repmat(DataM.Gen.DFIG.P_mec(:,Et), [lenWn 1]);
		DataM.Gen.DFIG.PQnormIE = squeeze(DataM.Gen.DFIG.PQnorm(1,2,:)) * one;
		DataM.Gen.DFIG.PQnormIF = squeeze(DataM.Gen.DFIG.PQnorm(1,3,:)) * one;
		DataM.Gen.DFIG.rIE = squeeze(DataM.Gen.DFIG.r(1,2,:)) * one;
		DataM.Gen.DFIG.rIF = squeeze(DataM.Gen.DFIG.r(1,3,:)) * one;
		DataM.Gen.DFIG.rOR = squeeze(DataM.Gen.DFIG.r(4,5,:)) * one;
		DataM.Gen.DFIG.sTopF = DataM.Gen.DFIG.sTop(3,:)' * one;
		DataM.Gen.DFIG.sTopR = DataM.Gen.DFIG.sTop(5,:)' * one;
		DataM.Gen.DFIG.xIE = squeeze(DataM.Gen.DFIG.x(1,2,:)) * one;
		DataM.Gen.DFIG.xIF = squeeze(DataM.Gen.DFIG.x(1,3,:)) * one;
		DataM.Gen.DFIG.xOR = squeeze(DataM.Gen.DFIG.x(4,5,:)) * one;
		DataM.Gen.DFIG.xiTopF = DataM.Gen.DFIG.xiTop(3,:)' * one;
		DataM.Gen.DFIG.xiTopR = DataM.Gen.DFIG.xiTop(5,:)' * one;
		DataM.Gen.DFIG.uLowE = DataM.Gen.DFIG.uLow(2,:)' * one;
		DataM.Gen.DFIG.uLowF = DataM.Gen.DFIG.uLow(3,:)' * one;
		DataM.Gen.DFIG.uTopE = DataM.Gen.DFIG.uTop(2,:)' * one;
		DataM.Gen.DFIG.uTopF = DataM.Gen.DFIG.uTop(3,:)' * one;
	end

	DataM.Gen.Pv.cr = DataM.Gen.Pv.cr * one;
	DataM.Gen.Pv.cv = DataM.Gen.Pv.cv * one;
	DataM.Gen.Pv.I = DataM.Gen.Pv.I * one;
	DataM.Gen.Pv.pPvg = DataM.Gen.Pv.pPvg * one;
	DataM.Gen.Pv.sTop = DataM.Gen.Pv.sTop * one;
	DataM.Gen.Pv.xiTop = DataM.Gen.Pv.xiTop * one;

	DataM.Gen.Tras.pgLow = DataM.Gen.Tras.pgLow(:,(1:lenEt));
	DataM.Gen.Tras.pgTop = DataM.Gen.Tras.pgTop(:,(1:lenEt));
	DataM.Gen.Tras.qgLow = DataM.Gen.Tras.qgLow(:,(1:lenEt));
	DataM.Gen.Tras.qgTop = DataM.Gen.Tras.qgTop(:,(1:lenEt));

	
	
	
	DataM.St.AC.beta = DataM.St.AC.beta * one;
	DataM.St.AC.epsilon = DataM.St.AC.epsilon * one;
	DataM.St.AC.eta = DataM.St.AC.eta * one;
	DataM.St.AC.tempLow = DataM.St.AC.tempLow * one;
	DataM.St.AC.tempPref = DataM.St.AC.tempPref * one;
	DataM.St.AC.tempTop = DataM.St.AC.tempTop * one;

	DataM.St.Bat.beta = DataM.St.Bat.beta * one;
	DataM.St.Bat.cr = DataM.St.Bat.cr * one;
	DataM.St.Bat.cv = DataM.St.Bat.cv * one;
	DataM.St.Bat.ELow = DataM.St.Bat.ELow * one;
	DataM.St.Bat.ETop = DataM.St.Bat.ETop * one;
	DataM.St.Bat.epsilon = DataM.St.Bat.epsilon * one;
	DataM.St.Bat.eta = DataM.St.Bat.eta * one;
	DataM.St.Bat.I = DataM.St.Bat.I * one;
	DataM.St.Bat.m1 = DataM.St.Bat.m1 * one;
	DataM.St.Bat.m2 = DataM.St.Bat.m2 * one;
	DataM.St.Bat.m3 = DataM.St.Bat.m3 * one;
	DataM.St.Bat.pgLow = DataM.St.Bat.pgLow * one;
	DataM.St.Bat.pgTop = DataM.St.Bat.pgTop * one;
	DataM.St.Bat.sTop = DataM.St.Bat.sTop * one;
	DataM.St.Bat.xiTop = DataM.St.Bat.xiTop * one;
	DataM.St.Bat.wOm = DataM.St.Bat.wOm * one;
	DataM.St.Bat.wU = DataM.St.Bat.wU * one;

	DataM.temp = DataM.temp * one;

	DataM.Util.betaE = DataM.Util.betaE(:,(1:lenEt));
	DataM.Util.betaT = DataM.Util.betaT * one;
	DataM.Util.pzCnLow = DataM.Util.pzCnLow(:,(1:lenEt),:);
	DataM.Util.pzCnLowE = DataM.Util.pzCnLowE * one;
	DataM.Util.pzCnPref = DataM.Util.pzCnPref(:,(1:lenEt),:);
	DataM.Util.pzCnPrefE = DataM.Util.pzCnPrefE * one;
	DataM.Util.pzCnTop = DataM.Util.pzCnTop(:,(1:lenEt),:);
	DataM.Util.pzCnTopE = DataM.Util.pzCnTopE * one;
	DataM.Util.qzCnLowE = DataM.Util.qzCnLowE * one;
	DataM.Util.qzCnTopE = DataM.Util.qzCnTopE * one;

end