function [DataM] = reshapeData(Data, Config)

	Buses = length(Data.Red.Branch.T);

	Data.Gen.DFIG.r	 = full(	Data.Gen.DFIG.r	);
	Data.Gen.DFIG.x	 = full(	Data.Gen.DFIG.x	);
	Data.Gen.DFIG.Tg	 = full(	Data.Gen.DFIG.Tg	);
	Data.Gen.DFIG.lTop	 = full(	Data.Gen.DFIG.lTop	);
	Data.Gen.DFIG.sTop	 = full(	Data.Gen.DFIG.sTop	);
	Data.temp	 = 	Data.temp(1:Config.Etapas);
	Data.Util.pzCnLow	 = 	Data.Util.pzCnLow(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1),:);
	Data.Util.pzCnPref	 = 	Data.Util.pzCnPref(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1),:);
	Data.Util.pzCnTop	 = 	Data.Util.pzCnTop(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1),:);
    Data.Util.betaE      =  Data.Util.betaE(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
    
% 	Data.Util.pzCnPrefE	 = 	Data.Util.pzCnPrefE(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
% 	Data.Util.pzCnLowE	 = 	Data.Util.pzCnLowE(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
% 	Data.Util.pzCnTopE	 = 	Data.Util.pzCnTopE(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));

	Data.Cost.piPTras	 = 	Data.Cost.piPTras(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
	Data.Cost.piQmtras	 = 	Data.Cost.piQmtras(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
	Data.Cost.piQMtras	 = 	Data.Cost.piQMtras(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
	Data.Cost.m	 = 	Data.Cost.m(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
	Data.Gen.DFIG.n_	 = 	Data.Gen.DFIG.n_(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
	Data.Gen.DFIG.P_mec	 = 	Data.Gen.DFIG.P_mec(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
	Data.Red.Bus.pCLow	 = 	Data.Red.Bus.pCLow(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));
	Data.Red.Bus.qCLow	 = 	Data.Red.Bus.qCLow(:,(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1));

	Tup = triu(Data.Red.Branch.T);
% 	Tdown = Data.Red.Branch.T .* (1 - Tup);

	
	DataM = Data;

	% 4D
	DataM.Red.Bus.alpha	 = 	replicateMat4_4(	Data.Red.Bus.alpha	,Config.Etapas);
	DataM.St.AC.a	 = 	replicateMat4_4(	Data.St.AC.a	,Config.Etapas);
	DataM.St.AC.beta	 = 	replicateMat4_4(	Data.St.AC.beta	,Config.Etapas);
	DataM.St.AC.epsilon	 = 	replicateMat4_4(	Data.St.AC.epsilon	,Config.Etapas);
	DataM.St.AC.eta	 = 	replicateMat4_4(	Data.St.AC.eta	,Config.Etapas);
	DataM.St.AC.tempLow	 = 	replicateMat4_4(	Data.St.AC.tempLow	,Config.Etapas);
	DataM.St.AC.tempPref	 = 	replicateMat4_4(	Data.St.AC.tempPref	,Config.Etapas);
	DataM.St.AC.tempTop	 = 	replicateMat4_4(	Data.St.AC.tempTop	,Config.Etapas);
	DataM.Util.aT	 = 	replicateMat4_3_4(	Data.Util.aT	,Config.Etapas,2);
	DataM.Util.betaT	 = 	replicateMat4_3_4(	Data.Util.betaT	,Config.Etapas,2);
	DataM.Util.pzCnLow	 = 	replicateMat4_2(	Data.Util.pzCnLow	,1);
	DataM.Util.pzCnTop	 = 	replicateMat4_2(	Data.Util.pzCnTop	,1);
	DataM.Util.pzCnPref	 = 	replicateMat4_2(	Data.Util.pzCnPref	,1);

	% DistrInfoM.Cl.v	 = 	replicateMat4DApp(	DistrInfo.Cl.v	,2);
	% DistrInfoM.lambdaT	 = 	replicateMat4DApp(	DistrInfo.lambdaT	,2);
	% DistrInfoM.muT	 = 	replicateMat4DApp(	DistrInfo.muT	,2);
    DataM.Gen.DFIG.r	= replicateMat4_3(	Data.Gen.DFIG.r	, Config.Etapas);
    DataM.Gen.DFIG.x	= replicateMat4_3(	Data.Gen.DFIG.x	, Config.Etapas);
    DataM.Gen.DFIG.lTop	= replicateMat4_3(	Data.Gen.DFIG.lTop	, Config.Etapas);
    DataM.Gen.DFIG.PQnorm	= replicateMat4_3(	Data.Gen.DFIG.PQnorm	, Config.Etapas);

    DataM.Gen.DFIG.sTop	= replicateMat4_3(	Data.Gen.DFIG.sTop	, Config.Etapas);

    DataM.Gen.DFIG.uLow	= replicateMat4_3(	Data.Gen.DFIG.uLow	, Config.Etapas);
    DataM.Gen.DFIG.uTop	= replicateMat4_3(	Data.Gen.DFIG.uTop	, Config.Etapas);
    DataM.Gen.DFIG.cv	= replicateMat4_3(	Data.Gen.DFIG.cv	, Config.Etapas);
    DataM.Gen.DFIG.cr	= replicateMat4_3(	Data.Gen.DFIG.cr	, Config.Etapas);

	% 3D
	DataM.temp	 = 	replicateMat3_1(	Data.temp'	,Buses);
	DataM.Cost.rhomqPv	 = 	replicateMat3_3(	Data.Cost.rhomqPv	,Config.Etapas);
	DataM.Cost.rhoMqPv	 = 	replicateMat3_3(	Data.Cost.rhoMqPv	,Config.Etapas);
	DataM.Gen.Pv.pPvg	 = 	replicateMat3_3(	Data.Gen.Pv.pPvg	,Config.Etapas);
	DataM.Gen.Pv.cv	 = 	replicateMat3_3(	Data.Gen.Pv.cv	,Config.Etapas);
	DataM.Gen.Pv.cr	 = 	replicateMat3_3(	Data.Gen.Pv.cr	,Config.Etapas);
	DataM.Gen.Pv.sTop	 = 	replicateMat3_3(	Data.Gen.Pv.sTop	,Config.Etapas);
	DataM.Cost.rhopPv	 = 	replicateMat3_3(	Data.Cost.rhopPv	,Config.Etapas);
	DataM.St.Bat.beta	 = 	replicateMat3_3(	Data.St.Bat.beta 	,Config.Etapas);
	DataM.St.Bat.cr	 = 	replicateMat3_3(	Data.St.Bat.cr 	,Config.Etapas);
	DataM.St.Bat.cv	 = 	replicateMat3_3(	Data.St.Bat.cv 	,Config.Etapas);
	DataM.St.Bat.EIni	 = 	replicateMat3_3(	Data.St.Bat.EIni 	,Config.Etapas);
	DataM.St.Bat.ELow	 = 	replicateMat3_3(	Data.St.Bat.ELow 	,Config.Etapas);
	DataM.St.Bat.epsilon	 = 	replicateMat3_3(	Data.St.Bat.epsilon 	,Config.Etapas);
	DataM.St.Bat.eta	 = 	replicateMat3_3(	Data.St.Bat.eta 	,Config.Etapas);
	DataM.St.Bat.ETop	 = 	replicateMat3_3(	Data.St.Bat.ETop 	,Config.Etapas);
	DataM.St.Bat.m1	 = 	replicateMat3_3(	Data.St.Bat.m1 	,Config.Etapas);
	DataM.St.Bat.m2	 = 	replicateMat3_3(	Data.St.Bat.m2 	,Config.Etapas);
	DataM.St.Bat.m3	 = 	replicateMat3_3(	Data.St.Bat.m3 	,Config.Etapas);
	DataM.St.Bat.pgLow	 = 	replicateMat3_3(	Data.St.Bat.pgLow 	,Config.Etapas);
	DataM.St.Bat.pgTop	 = 	replicateMat3_3(	Data.St.Bat.pgTop 	,Config.Etapas);
	DataM.St.Bat.sTop	 = 	replicateMat3_3(	Data.St.Bat.sTop 	,Config.Etapas);
	DataM.St.Bat.wOm	 = 	replicateMat3_3(	Data.St.Bat.wOm 	,Config.Etapas);
	DataM.St.Bat.wU	 = 	replicateMat3_3(	Data.St.Bat.wU 	,Config.Etapas);
	DataM.Cost.piPTras	 = 	replicateMat3_2(	Data.Cost.piPTras	,1);
	DataM.Cost.piQmtras	 = 	replicateMat3_2(	Data.Cost.piQmtras	,1);
	DataM.Cost.piQMtras	 = 	replicateMat3_2(	Data.Cost.piQMtras	,1);
	DataM.Red.Bus.P0Low	 = 	ones(Buses,1,Config.Etapas)*	Data.Red.Bus.P0Low	;
	DataM.Red.Bus.P0Top	 = 	ones(Buses,1,Config.Etapas)*	Data.Red.Bus.P0Top	;
	DataM.Red.Bus.Q0Low	 = 	ones(Buses,1,Config.Etapas)*	Data.Red.Bus.Q0Low	;
	DataM.Red.Bus.Q0Top	 = 	ones(Buses,1,Config.Etapas)*	Data.Red.Bus.Q0Top	;
	DataM.Cost.m	 = 	replicateMat3_1(	Data.Cost.m	,Buses);
	DataM.Gen.DFIG.n_	 = 	replicateMat3_1(	Data.Gen.DFIG.n_	,Buses);
	DataM.Gen.DFIG.P_mec	 = 	replicateMat3_1(	Data.Gen.DFIG.P_mec	,Buses);
	DataM.Gen.DFIG.Tg	 = 	replicateMat3_3(	Data.Gen.DFIG.Tg	, Config.Etapas);
    DataM.Gen.DFIG.I	= replicateMat3_3(	Data.Gen.DFIG.I	, Config.Etapas);
	DataM.Red.Bus.CapLow	 = 	replicateMat3_3(	Data.Red.Bus.CapLow	,Config.Etapas);
	DataM.Red.Bus.CapTop	 = 	replicateMat3_3(	Data.Red.Bus.CapTop	,Config.Etapas);
	DataM.Red.Bus.Ncp	 = 	replicateMat3_3(	Data.Red.Bus.Ncp	,Config.Etapas);
	DataM.Red.Bus.Ntr	 = 	replicateMat3_3(	Data.Red.Bus.Ntr	,Config.Etapas);
	DataM.Cost.cdv	 = 	replicateMat3_3(	Data.Cost.cdv	,Config.Etapas);
	DataM.Red.Bus.TapLow	 = 	replicateMat3_3(	Data.Red.Bus.TapLow	,Config.Etapas);
	DataM.Red.Bus.TapTop	 = 	replicateMat3_3(	Data.Red.Bus.TapTop	,Config.Etapas);
	DataM.Red.Bus.uLow	 = 	replicateMat3_3(	Data.Red.Bus.uLow	,Config.Etapas);
	DataM.Red.Bus.uTop	 = 	replicateMat3_3(	Data.Red.Bus.uTop	,Config.Etapas);
	DataM.Red.Branch.lTop	 = 	replicateMat3_3(	Data.Red.Branch.lTop 	,Config.Etapas);
	DataM.Red.Branch.r	 = 	replicateMat3_3(	Data.Red.Branch.r	,Config.Etapas);
	DataM.Red.Branch.T	 = 	replicateMat3_3(	Data.Red.Branch.T	,Config.Etapas);
	DataM.Red.Branch.Tup	 = 	replicateMat3_3(	Tup	,Config.Etapas);
% 	DataM.Red.Branch.Tdown	 = 	replicateMat3_3(	Tdown	,Config.Etapas);
	DataM.Red.Branch.x	 = 	replicateMat3_3(	Data.Red.Branch.x	,Config.Etapas);
	DataM.Red.Branch.yLow	 = 	replicateMat3_3(	Data.Red.Branch.yLow	,Config.Etapas);
	DataM.Red.Branch.yTop	 = 	replicateMat3_3(	Data.Red.Branch.yTop	,Config.Etapas);
	DataM.Red.Bus.pCLow	 = 	replicateMat3_2(	Data.Red.Bus.pCLow	,1);
	DataM.Red.Bus.qCLow	 = 	replicateMat3_2(	Data.Red.Bus.qCLow	,1);
	
	DataM.Util.aE	 = 	replicateMat3_3(	Data.Util.aE	,Config.Etapas);
	DataM.Util.betaE	 = 	replicateMat3_2(	Data.Util.betaE	,1);

	DataM.Util.pzCnPrefE	 = 	replicateMat3_3(	Data.Util.pzCnPrefE	,Config.Etapas);
	DataM.Util.pzCnLowE	 = 	replicateMat3_3(	Data.Util.pzCnLowE	,Config.Etapas);
	DataM.Util.pzCnTopE	 = 	replicateMat3_3(	Data.Util.pzCnTopE	,Config.Etapas);

	DataM.St.Bat.I	 = 	replicateMat3_3(	Data.St.Bat.I	,Config.Etapas);
	DataM.Gen.Pv.I	 = 	replicateMat3_3(	Data.Gen.Pv.I	,Config.Etapas);

	DataM.ClNI.pC	 = 	replicateMat3_3(	Data.ClNI.pC	,Config.Etapas);
	DataM.ClNI.qC	 = 	replicateMat3_3(	Data.ClNI.qC	,Config.Etapas);


    DataM.Cost.rhopWi	= replicateMat3_3(	Data.Cost.rhopWi	, Config.Etapas);
    DataM.Cost.rhomqWi	= replicateMat3_3(	Data.Cost.rhomqWi	, Config.Etapas);
    DataM.Cost.rhoMqWi	= replicateMat3_3(	Data.Cost.rhoMqWi	, Config.Etapas);
    
    

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
