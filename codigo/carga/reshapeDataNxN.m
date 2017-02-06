function [DataNxN] = reshapeDataNxN(Data, Config)

	Buses = length(Data.Red.Branch.T);
	iniEt = Config.iniEtapa;
	finEt = iniEt+Config.Etapas-1;
	Et = (iniEt:finEt);
	
	DataNxN = Data;

	Tup = triu(DataNxN.Red.Branch.T);

	
	DataNxN.Red.Branch.cY	 = 	replicateMat3_3(	DataNxN.Red.Branch.cY	,Config.Etapas);
	DataNxN.Red.Branch.lTop	 = 	replicateMat3_3(	DataNxN.Red.Branch.lTop 	,Config.Etapas);
	DataNxN.Red.Branch.r	 = 	replicateMat3_3(	DataNxN.Red.Branch.r	,Config.Etapas);
	DataNxN.Red.Branch.T	 = 	replicateMat3_3(	DataNxN.Red.Branch.T	,Config.Etapas);
	DataNxN.Red.Branch.Tup	 = 	replicateMat3_3(	Tup	,Config.Etapas);
	DataNxN.Red.Branch.x	 = 	replicateMat3_3(	DataNxN.Red.Branch.x	,Config.Etapas);
	DataNxN.Red.Branch.yLow	 = 	replicateMat3_3(	DataNxN.Red.Branch.yLow	,Config.Etapas);
	DataNxN.Red.Branch.yTop	 = 	replicateMat3_3(	DataNxN.Red.Branch.yTop	,Config.Etapas);
	
	% 4D
	DataNxN.Red.Bus.alpha	 = 	replicateMat4_4(	DataNxN.Red.Bus.alpha	,Config.Etapas);
	% 3D
	DataNxN.Red.Bus.CapLow	 = 	replicateMat3_3(	DataNxN.Red.Bus.CapLow	,Config.Etapas);
	DataNxN.Red.Bus.CapTop	 = 	replicateMat3_3(	DataNxN.Red.Bus.CapTop	,Config.Etapas);
	DataNxN.Red.Bus.indCap	 = 	replicateMat3_3(	DataNxN.Red.Bus.indCap	,Config.Etapas);
	DataNxN.Red.Bus.indTap	 = 	replicateMat3_3(	DataNxN.Red.Bus.indTap	,Config.Etapas);
	DataNxN.Red.Bus.Ncp	 = 	replicateMat3_3(	DataNxN.Red.Bus.Ncp	,Config.Etapas);
	DataNxN.Red.Bus.Ntr	 = 	replicateMat3_3(	DataNxN.Red.Bus.Ntr	,Config.Etapas);
	DataNxN.Red.Bus.pCLow	 = 	replicateMat3_2(	DataNxN.Red.Bus.pCLow(:,Et)	,1);
	DataNxN.Red.Bus.qCLow	 = 	replicateMat3_2(	DataNxN.Red.Bus.qCLow(:,Et)	,1);
	DataNxN.Red.Bus.TapLow	 = 	replicateMat3_3(	DataNxN.Red.Bus.TapLow	,Config.Etapas);
	DataNxN.Red.Bus.TapTop	 = 	replicateMat3_3(	DataNxN.Red.Bus.TapTop	,Config.Etapas);
	DataNxN.Red.Bus.uLow	 = 	replicateMat3_3(	DataNxN.Red.Bus.uLow	,Config.Etapas);
	DataNxN.Red.Bus.uTop	 = 	replicateMat3_3(	DataNxN.Red.Bus.uTop	,Config.Etapas);

	
	DataNxN.ClNI.pC	 = 	replicateMat3_3(	DataNxN.ClNI.pC	,Config.Etapas);
	DataNxN.ClNI.qC	 = 	replicateMat3_3(	DataNxN.ClNI.qC	,Config.Etapas);


	
	
	
	DataNxN.Cost.cdv	 = 	replicateMat3_3(	DataNxN.Cost.cdv	,Config.Etapas);
    DataNxN.Cost.m	 = 	replicateMat3_1(	DataNxN.Cost.m(:,Et)	,Buses);
	DataNxN.Cost.piPTras	 = 	replicateMat3_2(	DataNxN.Cost.piPTras(:,Et)	,1);
	DataNxN.Cost.piQmtras	 = 	replicateMat3_2(	DataNxN.Cost.piQmtras(:,Et)	,1);
	DataNxN.Cost.piQMtras	 = 	replicateMat3_2(	DataNxN.Cost.piQMtras(:,Et)	,1);
	DataNxN.Cost.rhopPv	 = 	replicateMat3_3(	DataNxN.Cost.rhopPv	,Config.Etapas);
	DataNxN.Cost.rhomqPv	 = 	replicateMat3_3(	DataNxN.Cost.rhomqPv	,Config.Etapas);
	DataNxN.Cost.rhoMqPv	 = 	replicateMat3_3(	DataNxN.Cost.rhoMqPv	,Config.Etapas);
    DataNxN.Cost.rhopWi	= replicateMat3_3(	DataNxN.Cost.rhopWi	, Config.Etapas);
    DataNxN.Cost.rhomqWi	= replicateMat3_3(	DataNxN.Cost.rhomqWi	, Config.Etapas);
    DataNxN.Cost.rhoMqWi	= replicateMat3_3(	DataNxN.Cost.rhoMqWi	, Config.Etapas);

	
    DataNxN.Gen.DFIG.cr	= replicateMat4_3(	DataNxN.Gen.DFIG.cr	, Config.Etapas);
    DataNxN.Gen.DFIG.cv	= replicateMat4_3(	DataNxN.Gen.DFIG.cv	, Config.Etapas);
    DataNxN.Gen.DFIG.I	= replicateMat3_3(	DataNxN.Gen.DFIG.I	, Config.Etapas);
    DataNxN.Gen.DFIG.lTop	= replicateMat4_3(full(	DataNxN.Gen.DFIG.lTop	), Config.Etapas);
	DataNxN.Gen.DFIG.n_	 = 	replicateMat3_1(	DataNxN.Gen.DFIG.n_(:,Et)	,Buses);
	DataNxN.Gen.DFIG.P_mec	 = 	replicateMat3_1(	DataNxN.Gen.DFIG.P_mec(:,Et)	,Buses);
    DataNxN.Gen.DFIG.PQnorm	= replicateMat4_3(	DataNxN.Gen.DFIG.PQnorm	, Config.Etapas);
    DataNxN.Gen.DFIG.r	= replicateMat4_3(full(	DataNxN.Gen.DFIG.r	), Config.Etapas);
    DataNxN.Gen.DFIG.sTop	= replicateMat4_3(full(	DataNxN.Gen.DFIG.sTop	), Config.Etapas);
	DataNxN.Gen.DFIG.Tg	 = 	replicateMat3_3(full(	DataNxN.Gen.DFIG.Tg	), Config.Etapas);
    DataNxN.Gen.DFIG.x	= replicateMat4_3(full(	DataNxN.Gen.DFIG.x	), Config.Etapas);
    DataNxN.Gen.DFIG.xiTop	= replicateMat4_3(full(	DataNxN.Gen.DFIG.xiTop	), Config.Etapas);
    DataNxN.Gen.DFIG.uLow	= replicateMat4_3(	DataNxN.Gen.DFIG.uLow	, Config.Etapas);
    DataNxN.Gen.DFIG.uTop	= replicateMat4_3(	DataNxN.Gen.DFIG.uTop	, Config.Etapas);

	DataNxN.Gen.Pv.cv	 = 	replicateMat3_3(	DataNxN.Gen.Pv.cv	,Config.Etapas);
	DataNxN.Gen.Pv.cr	 = 	replicateMat3_3(	DataNxN.Gen.Pv.cr	,Config.Etapas);
	DataNxN.Gen.Pv.I	 = 	replicateMat3_3(	DataNxN.Gen.Pv.I	,Config.Etapas);
	DataNxN.Gen.Pv.pPvg	 = 	replicateMat3_3(	DataNxN.Gen.Pv.pPvg	,Config.Etapas);
	DataNxN.Gen.Pv.sTop	 = 	replicateMat3_3(	DataNxN.Gen.Pv.sTop	,Config.Etapas);
	DataNxN.Gen.Pv.xiTop	 = 	replicateMat3_3(	DataNxN.Gen.Pv.xiTop	,Config.Etapas);

    DataNxN.Gen.Tras.pgLow = replicateMat3_2(	DataNxN.Gen.Tras.pgLow(:,Et)	,1);
	DataNxN.Gen.Tras.pgTop = replicateMat3_2(	DataNxN.Gen.Tras.pgTop(:,Et)	,1);
	DataNxN.Gen.Tras.qgLow = replicateMat3_2(	DataNxN.Gen.Tras.qgLow(:,Et)	,1);
	DataNxN.Gen.Tras.qgTop = replicateMat3_2(	DataNxN.Gen.Tras.qgTop(:,Et)	,1);


	% DataNxN.St.AC.a	 = 	replicateMat4_4(	DataNxN.St.AC.a	,Config.Etapas);
	DataNxN.St.AC.beta	 = 	replicateMat4_4(	DataNxN.St.AC.beta	,Config.Etapas);
	DataNxN.St.AC.epsilon	 = 	replicateMat4_4(	DataNxN.St.AC.epsilon	,Config.Etapas);
	DataNxN.St.AC.eta	 = 	replicateMat4_4(	DataNxN.St.AC.eta	,Config.Etapas);
	DataNxN.St.AC.tempLow	 = 	replicateMat4_4(	DataNxN.St.AC.tempLow	,Config.Etapas);
	DataNxN.St.AC.tempPref	 = 	replicateMat4_4(	DataNxN.St.AC.tempPref	,Config.Etapas);
	DataNxN.St.AC.tempTop	 = 	replicateMat4_4(	DataNxN.St.AC.tempTop	,Config.Etapas);
	
	DataNxN.St.Bat.beta	 = 	replicateMat3_3(	DataNxN.St.Bat.beta 	,Config.Etapas);
	DataNxN.St.Bat.cr	 = 	replicateMat3_3(	DataNxN.St.Bat.cr 	,Config.Etapas);
	DataNxN.St.Bat.cv	 = 	replicateMat3_3(	DataNxN.St.Bat.cv 	,Config.Etapas);
	DataNxN.St.Bat.EIni	 = 	replicateMat3_3(	DataNxN.St.Bat.EIni 	,Config.Etapas);
	DataNxN.St.Bat.ELow	 = 	replicateMat3_3(	DataNxN.St.Bat.ELow 	,Config.Etapas);
	DataNxN.St.Bat.epsilon	 = 	replicateMat3_3(	DataNxN.St.Bat.epsilon 	,Config.Etapas);
	DataNxN.St.Bat.eta	 = 	replicateMat3_3(	DataNxN.St.Bat.eta 	,Config.Etapas);
	DataNxN.St.Bat.ETop	 = 	replicateMat3_3(	DataNxN.St.Bat.ETop 	,Config.Etapas);
	DataNxN.St.Bat.I	 = 	replicateMat3_3(	DataNxN.St.Bat.I	,Config.Etapas);
	DataNxN.St.Bat.m1	 = 	replicateMat3_3(	DataNxN.St.Bat.m1 	,Config.Etapas);
	DataNxN.St.Bat.m2	 = 	replicateMat3_3(	DataNxN.St.Bat.m2 	,Config.Etapas);
	DataNxN.St.Bat.m3	 = 	replicateMat3_3(	DataNxN.St.Bat.m3 	,Config.Etapas);
	DataNxN.St.Bat.pgLow	 = 	replicateMat3_3(	DataNxN.St.Bat.pgLow 	,Config.Etapas);
	DataNxN.St.Bat.pgTop	 = 	replicateMat3_3(	DataNxN.St.Bat.pgTop 	,Config.Etapas);
	DataNxN.St.Bat.sTop	 = 	replicateMat3_3(	DataNxN.St.Bat.sTop 	,Config.Etapas);
	DataNxN.St.Bat.xiTop	 = 	replicateMat3_3(	DataNxN.St.Bat.xiTop 	,Config.Etapas);
	DataNxN.St.Bat.wOm	 = 	replicateMat3_3(	DataNxN.St.Bat.wOm 	,Config.Etapas);
	DataNxN.St.Bat.wU	 = 	replicateMat3_3(	DataNxN.St.Bat.wU 	,Config.Etapas);

	DataNxN.temp	 = 	replicateMat3_1(	DataNxN.temp(Et)'	,Buses);

	% DataNxN.Util.aT	 = 	replicateMat4_3_4(	DataNxN.Util.aT	,Config.Etapas,2);
	% DataNxN.Util.aE	 = 	replicateMat3_3(	DataNxN.Util.aE	,Config.Etapas);
	DataNxN.Util.betaE	 = 	replicateMat3_2(	DataNxN.Util.betaE(:,Et)	,1);
	DataNxN.Util.betaT	 = 	replicateMat4_3_4(	DataNxN.Util.betaT	,Config.Etapas,2);
	DataNxN.Util.pzCnLow	 = 	replicateMat4_2(	DataNxN.Util.pzCnLow(:,Et,:)	,1);
	DataNxN.Util.pzCnLowE	 = 	replicateMat3_3(	DataNxN.Util.pzCnLowE	,Config.Etapas);
	DataNxN.Util.pzCnPref	 = 	replicateMat4_2(	DataNxN.Util.pzCnPref(:,Et,:)	,1);
	DataNxN.Util.pzCnPrefE	 = 	replicateMat3_3(	DataNxN.Util.pzCnPrefE	,Config.Etapas);
	DataNxN.Util.pzCnTop	 = 	replicateMat4_2(	DataNxN.Util.pzCnTop(:,Et,:)	,1);
	DataNxN.Util.pzCnTopE	 = 	replicateMat3_3(	DataNxN.Util.pzCnTopE	,Config.Etapas);
	DataNxN.Util.qzCnPrefE	 = 	replicateMat3_3(	DataNxN.Util.qzCnPrefE	,Config.Etapas);
	DataNxN.Util.qzCnLowE	 = 	replicateMat3_3(	DataNxN.Util.qzCnLowE	,Config.Etapas);
	DataNxN.Util.qzCnTopE	 = 	replicateMat3_3(	DataNxN.Util.qzCnTopE	,Config.Etapas);



    
    

end

function [A] = replicateMat4_4(V,D4)
	A = zeros(size(V,1),1,D4,size(V,2));
	for d = 1:D4
		A(:,1,d,:) = V;
	end
end

function [A] = replicateMat4_3_4(V,D3,D4)
	A = zeros(size(V,1),1,D3,D4);
	for d3 = 1:D3
		for d4 = 1:D4
			A(:,1,d3,d4) = V;
		end
	end
end

function [A] = replicateMat4_3(V,D3)
	A = zeros(size(V,1),size(V,2),D3,size(V,3));
	for d3 = 1:D3
			A(:,:,d3,:) = V;
	end
end

function [A] = replicateMat4_2(V,D2)
	A = zeros(size(V,1),D2,size(V,2),size(V,3));
	for d = 1:D2
		A(:,d,:,:) = V;
	end
end

%-----------------------------------------------------

function [A] = replicateMat3_3(V,D3)
	A = zeros(size(V,1),size(V,2),D3);
	for d = 1:D3
		A(:,:,d) = V;
	end
end

function [A] = replicateMat3_2(V,D2)
	A = zeros(size(V,1),D2,size(V,2));
	for d = 1:D2
		A(:,d,:) = V;
	end
end

function [A] = replicateMat3_1(V,D1)
	A = zeros(D1,size(V,1),size(V,2));
	for d = 1:D1
		A(d,:,:) = V;
	end
end
