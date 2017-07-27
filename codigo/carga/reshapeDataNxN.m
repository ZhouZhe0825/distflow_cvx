function [DataNxN] = reshapeDataNxN(Data, Config)

	Buses = length(Data.Red.Branch.T);
	iniEt = Config.iniEtapa;
	finEt = iniEt+Config.Etapas-1;
	Et = (iniEt:finEt);
	
	DataNxN = Data;

	Tup = triu(DataNxN.Red.Branch.T);

	
	DataNxN.Red.Branch.lTop	 = 	replicateMat3_3(	DataNxN.Red.Branch.lTop 	,Config.Etapas);
	DataNxN.Red.Branch.r	 = 	replicateMat3_3(	DataNxN.Red.Branch.r	,Config.Etapas);
	DataNxN.Red.Branch.T	 = 	replicateMat3_3(	DataNxN.Red.Branch.T	,Config.Etapas);
	DataNxN.Red.Branch.Tup	 = 	replicateMat3_3(	Tup	,Config.Etapas);
	DataNxN.Red.Branch.x	 = 	replicateMat3_3(	DataNxN.Red.Branch.x	,Config.Etapas);
	DataNxN.Red.Branch.yLow	 = 	replicateMat3_3(	DataNxN.Red.Branch.yLow	,Config.Etapas);
	DataNxN.Red.Branch.yTop	 = 	replicateMat3_3(	DataNxN.Red.Branch.yTop	,Config.Etapas);
	DataNxN.Red.Branch.Tap	 = 	replicateMat3_3(	DataNxN.Red.Branch.Tap	,Config.Etapas);
	DataNxN.Red.Branch.NtrLow	 = 	replicateMat3_3(	DataNxN.Red.Branch.NtrLow	,Config.Etapas);
	DataNxN.Red.Branch.NtrTop	 = 	replicateMat3_3(	DataNxN.Red.Branch.NtrTop	,Config.Etapas);
	
	% 4D
	DataNxN.Red.Bus.alpha	 = 	replicateMat4_4(	DataNxN.Red.Bus.alpha	,Config.Etapas);
	% 3D
	DataNxN.Red.Bus.NcpLow	 = 	replicateMat3_3(	DataNxN.Red.Bus.NcpLow	,Config.Etapas);
	DataNxN.Red.Bus.NcpTop	 = 	replicateMat3_3(	DataNxN.Red.Bus.NcpTop	,Config.Etapas);
	DataNxN.Red.Bus.Cap	 = 	replicateMat3_3(	DataNxN.Red.Bus.Cap	,Config.Etapas);
	DataNxN.Red.Bus.pCLow	 = 	replicateMat3_2(	DataNxN.Red.Bus.pCLow(:,Et)	,1);
	DataNxN.Red.Bus.qCLow	 = 	replicateMat3_2(	DataNxN.Red.Bus.qCLow(:,Et)	,1);
	DataNxN.Red.Bus.uLow	 = 	replicateMat3_3(	DataNxN.Red.Bus.uLow	,Config.Etapas);
	DataNxN.Red.Bus.uTop	 = 	replicateMat3_3(	DataNxN.Red.Bus.uTop	,Config.Etapas);

	
	DataNxN.ClNI.pC	 = 	replicateMat3_3(	DataNxN.ClNI.pC	,Config.Etapas);
	DataNxN.ClNI.qC	 = 	replicateMat3_3(	DataNxN.ClNI.qC	,Config.Etapas);


	
	
	
	DataNxN.Cost.cdv	 = 	replicateMat3_2(	DataNxN.Cost.cdv(:,Et)	,1);
	DataNxN.Cost.cCap	 = 	replicateMat3_3(	DataNxN.Cost.cCap	,Config.Etapas);
	DataNxN.Cost.cTap	 = 	replicateMat3_3(	DataNxN.Cost.cTap	,Config.Etapas);
	DataNxN.Cost.cY	 = 	replicateMat3_3(	DataNxN.Cost.cY	,Config.Etapas);
	DataNxN.Cost.delta	 = 	replicateMat3_2(	DataNxN.Cost.delta(:,Et)	,1);
	DataNxN.Cost.m	 = 	replicateMat3_2(	DataNxN.Cost.m(:,Et)	,1);
	DataNxN.Cost.piPTras	 = 	replicateMat3_2(	DataNxN.Cost.piPTras(:,Et)	,1);
	DataNxN.Cost.piQmtras	 = 	replicateMat3_2(	DataNxN.Cost.piQmtras(:,Et)	,1);
	DataNxN.Cost.piQMtras	 = 	replicateMat3_2(	DataNxN.Cost.piQMtras(:,Et)	,1);
	DataNxN.Cost.rhopPv	 = 	replicateMat3_2(	DataNxN.Cost.rhopPv(:,Et)	,1);
	DataNxN.Cost.rhomqPv	 = 	replicateMat3_2(	DataNxN.Cost.rhomqPv(:,Et)	,1);
	DataNxN.Cost.rhoMqPv	 = 	replicateMat3_2(	DataNxN.Cost.rhoMqPv(:,Et)	,1);
	DataNxN.Cost.rhopWi	= replicateMat3_2(	DataNxN.Cost.rhopWi(:,Et)	, 1);
	DataNxN.Cost.rhomqWi	= replicateMat3_2(	DataNxN.Cost.rhomqWi(:,Et)	, 1);
	DataNxN.Cost.rhoMqWi	= replicateMat3_2(	DataNxN.Cost.rhoMqWi(:,Et)	, 1);

	
	% DataNxN.Gen.DFIG.cr	= replicateMat4_3(	DataNxN.Gen.DFIG.cr	, Config.Etapas);
	% DataNxN.Gen.DFIG.cv	= replicateMat4_3(	DataNxN.Gen.DFIG.cv	, Config.Etapas);
	% DataNxN.Gen.DFIG.I	= replicateMat3_3(	DataNxN.Gen.DFIG.I	, Config.Etapas);
	% DataNxN.Gen.DFIG.lTop	= replicateMat4_3(full(	DataNxN.Gen.DFIG.lTop	), Config.Etapas);
	DataNxN.Gen.DFIG.n_	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.n_(:,Et)	,1);
	DataNxN.Gen.DFIG.P_mec	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.P_mec(:,Et)	,1);
	% DataNxN.Gen.DFIG.PQnorm	= replicateMat4_3(	DataNxN.Gen.DFIG.PQnorm	, Config.Etapas);
	% DataNxN.Gen.DFIG.r	= replicateMat4_3(full(	DataNxN.Gen.DFIG.r	), Config.Etapas);
	% DataNxN.Gen.DFIG.sTop	= replicateMat4_3(full(	DataNxN.Gen.DFIG.sTop	), Config.Etapas);
% 	DataNxN.Gen.DFIG.Tg	 = 	replicateMat3_3(full(	DataNxN.Gen.DFIG.Tg	), Config.Etapas);
	% DataNxN.Gen.DFIG.x	= replicateMat4_3(full(	DataNxN.Gen.DFIG.x	), Config.Etapas);
	% DataNxN.Gen.DFIG.xiTop	= replicateMat4_3(full(	DataNxN.Gen.DFIG.xiTop	), Config.Etapas);
	% DataNxN.Gen.DFIG.uLow	= replicateMat4_3(	DataNxN.Gen.DFIG.uLow	, Config.Etapas);
	% DataNxN.Gen.DFIG.uTop	= replicateMat4_3(	DataNxN.Gen.DFIG.uTop	, Config.Etapas);


% 	cr	= replicateMat4_3(	DataNxN.Gen.DFIG.cr	, Config.Etapas);
% 	cv	= replicateMat4_3(	DataNxN.Gen.DFIG.cv	, Config.Etapas);
% 	I	= replicateMat3_3(	DataNxN.Gen.DFIG.I	, Config.Etapas);
% 	lTop	= replicateMat4_3(full(	DataNxN.Gen.DFIG.lTop	), Config.Etapas);
% 	PQnorm	= replicateMat4_3(	DataNxN.Gen.DFIG.PQnorm	, Config.Etapas);
% 	r	= replicateMat4_3(full(	DataNxN.Gen.DFIG.r	), Config.Etapas);
% 	sTop	= replicateMat4_3(full(	DataNxN.Gen.DFIG.sTop	), Config.Etapas);
% 	x	= replicateMat4_3(full(	DataNxN.Gen.DFIG.x	), Config.Etapas);
% 	xiTop	= replicateMat4_3(full(	DataNxN.Gen.DFIG.xiTop	), Config.Etapas);
% 	uLow	= replicateMat4_3(	DataNxN.Gen.DFIG.uLow	, Config.Etapas);
% 	uTop	= replicateMat4_3(	DataNxN.Gen.DFIG.uTop	, Config.Etapas);


    DataNxN.Gen.DFIG.C_plb	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.C_plb	,Config.Etapas);
    DataNxN.Gen.DFIG.G	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.G	,Config.Etapas);
    DataNxN.Gen.DFIG.N_er	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.N_er	,Config.Etapas);
    DataNxN.Gen.DFIG.Np	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.Np	,Config.Etapas);
    DataNxN.Gen.DFIG.Omega	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.Omega	,Config.Etapas);
    DataNxN.Gen.DFIG.PQnormIE	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.PQnormIE(:,Et)	,1);
    DataNxN.Gen.DFIG.PQnormIF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.PQnormIF(:,Et)	,1);
    DataNxN.Gen.DFIG.P_nMec	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.P_nMec	,Config.Etapas);
    DataNxN.Gen.DFIG.R_	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.R_	,Config.Etapas);
    DataNxN.Gen.DFIG.c1	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.c1	,Config.Etapas);
    DataNxN.Gen.DFIG.c2	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.c2	,Config.Etapas);
    DataNxN.Gen.DFIG.c3	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.c3	,Config.Etapas);
    DataNxN.Gen.DFIG.c4	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.c4	,Config.Etapas);
    DataNxN.Gen.DFIG.c5	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.c5	,Config.Etapas);
    DataNxN.Gen.DFIG.c6	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.c6	,Config.Etapas);
    DataNxN.Gen.DFIG.c7	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.c7	,Config.Etapas);
    DataNxN.Gen.DFIG.crF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.crF(:,Et)	,1);
    DataNxN.Gen.DFIG.crR	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.crR(:,Et)	,1);
    DataNxN.Gen.DFIG.cvF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.cvF(:,Et)	,1);
    DataNxN.Gen.DFIG.cvR	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.cvR(:,Et)	,1);
    DataNxN.Gen.DFIG.lTopIE	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.lTopIE(:,Et)	,1);
    DataNxN.Gen.DFIG.lTopIF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.lTopIF(:,Et)	,1);
    DataNxN.Gen.DFIG.lTopOR	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.lTopOR(:,Et)	,1);
    DataNxN.Gen.DFIG.lambda_opt	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.lambda_opt	,Config.Etapas);
    DataNxN.Gen.DFIG.rIE	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.rIE(:,Et)	,1);
    DataNxN.Gen.DFIG.rIF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.rIF(:,Et)	,1);
    DataNxN.Gen.DFIG.rOR	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.rOR(:,Et)	,1);
    DataNxN.Gen.DFIG.rho	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.rho	,Config.Etapas);
    DataNxN.Gen.DFIG.sTopF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.sTopF(:,Et)	,1);
    DataNxN.Gen.DFIG.sTopR	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.sTopR(:,Et)	,1);
    DataNxN.Gen.DFIG.uLowE	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.uLowE(:,Et)	,1);
    DataNxN.Gen.DFIG.uLowF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.uLowF(:,Et)	,1);
    DataNxN.Gen.DFIG.uTopE	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.uTopE(:,Et)	,1);
    DataNxN.Gen.DFIG.uTopF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.uTopF(:,Et)	,1);
    DataNxN.Gen.DFIG.vM	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.vM	,Config.Etapas);
    DataNxN.Gen.DFIG.vMM	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.vMM	,Config.Etapas);
    DataNxN.Gen.DFIG.vm	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.vm	,Config.Etapas);
    DataNxN.Gen.DFIG.vmm	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.vmm	,Config.Etapas);
    DataNxN.Gen.DFIG.vv	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.vv	,Config.Etapas);
    DataNxN.Gen.DFIG.w	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.w	,Config.Etapas);
    DataNxN.Gen.DFIG.ws	 = 	replicateMat3_3(	DataNxN.Gen.DFIG.ws	,Config.Etapas);
    DataNxN.Gen.DFIG.xIE	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.xIE(:,Et)	,1);
    DataNxN.Gen.DFIG.xIF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.xIF(:,Et)	,1);
    DataNxN.Gen.DFIG.xOR	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.xOR(:,Et)	,1);
    DataNxN.Gen.DFIG.xiTopF	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.xiTopF(:,Et)	,1);
    DataNxN.Gen.DFIG.xiTopR	 = 	replicateMat3_2(	DataNxN.Gen.DFIG.xiTopR(:,Et)	,1);


	DataNxN.Gen.Pv.cr	 = 	replicateMat3_3(	DataNxN.Gen.Pv.cr	,Config.Etapas);
	DataNxN.Gen.Pv.cv	 = 	replicateMat3_3(	DataNxN.Gen.Pv.cv	,Config.Etapas);
	DataNxN.Gen.Pv.pPvg	 = 	replicateMat3_2(	DataNxN.Gen.Pv.pPvg(:,Et)	,1);
	DataNxN.Gen.Pv.sTop	 = 	replicateMat3_3(	DataNxN.Gen.Pv.sTop	,Config.Etapas);
	DataNxN.Gen.Pv.xiTop	 = 	replicateMat3_3(	DataNxN.Gen.Pv.xiTop	,Config.Etapas);

    DataNxN.Gen.Tras.pgLow = replicateMat3_3(	DataNxN.Gen.Tras.pgLow	,Config.Etapas);
	DataNxN.Gen.Tras.pgTop = replicateMat3_3(	DataNxN.Gen.Tras.pgTop	,Config.Etapas);
	DataNxN.Gen.Tras.qgLow = replicateMat3_3(	DataNxN.Gen.Tras.qgLow	,Config.Etapas);
	DataNxN.Gen.Tras.qgTop = replicateMat3_3(	DataNxN.Gen.Tras.qgTop	,Config.Etapas);


	% DataNxN.St.AC.a	 = 	replicateMat4_4(	DataNxN.St.AC.a	,Config.Etapas);
	DataNxN.St.AC.beta	 = 	replicateMat3_2(	DataNxN.St.AC.beta(:,Et)	,1);
	DataNxN.St.AC.epsilon	 = 	replicateMat3_2(	DataNxN.St.AC.epsilon(:,Et)	,1);
	DataNxN.St.AC.eta	 = 	replicateMat3_2(	DataNxN.St.AC.eta(:,Et)	,1);
	DataNxN.St.AC.tempLow	 = 	replicateMat3_2(	DataNxN.St.AC.tempLow(:,Et)	,1);
	DataNxN.St.AC.tempPref	 = 	replicateMat3_2(	DataNxN.St.AC.tempPref(:,Et)	,1);
	DataNxN.St.AC.tempTop	 = 	replicateMat3_2(	DataNxN.St.AC.tempTop(:,Et)	,1);
	
	DataNxN.St.Bat.beta	 = 	replicateMat3_3(	DataNxN.St.Bat.beta 	,Config.Etapas);
	DataNxN.St.Bat.cr	 = 	replicateMat3_3(	DataNxN.St.Bat.cr 	,Config.Etapas);
	DataNxN.St.Bat.cv	 = 	replicateMat3_3(	DataNxN.St.Bat.cv 	,Config.Etapas);
	DataNxN.St.Bat.EIni	 = 	replicateMat3_3(	DataNxN.St.Bat.EIni 	,Config.Etapas);
	DataNxN.St.Bat.ELow	 = 	replicateMat3_3(	DataNxN.St.Bat.ELow 	,Config.Etapas);
	DataNxN.St.Bat.epsilon	 = 	replicateMat3_3(	DataNxN.St.Bat.epsilon 	,Config.Etapas);
	DataNxN.St.Bat.etaC	 = 	replicateMat3_3(	DataNxN.St.Bat.etaC 	,Config.Etapas);
	DataNxN.St.Bat.etaD	 = 	replicateMat3_3(	DataNxN.St.Bat.etaD 	,Config.Etapas);
	DataNxN.St.Bat.ETop	 = 	replicateMat3_3(	DataNxN.St.Bat.ETop 	,Config.Etapas);
	DataNxN.St.Bat.EPref	 = 	replicateMat3_3(	DataNxN.St.Bat.EPref 	,Config.Etapas);
	DataNxN.St.Bat.gama	 = 	replicateMat3_3(	DataNxN.St.Bat.gama 	,Config.Etapas);
	DataNxN.St.Bat.kapa	 = 	replicateMat3_3(	DataNxN.St.Bat.kapa 	,Config.Etapas);
	DataNxN.St.Bat.m1	 = 	replicateMat3_3(	DataNxN.St.Bat.m1 	,Config.Etapas);
	DataNxN.St.Bat.m2	 = 	replicateMat3_3(	DataNxN.St.Bat.m2 	,Config.Etapas);
	DataNxN.St.Bat.m3	 = 	replicateMat3_3(	DataNxN.St.Bat.m3 	,Config.Etapas);
	DataNxN.St.Bat.pgLow	 = 	replicateMat3_3(	DataNxN.St.Bat.pgLow 	,Config.Etapas);
	DataNxN.St.Bat.pgTop	 = 	replicateMat3_3(	DataNxN.St.Bat.pgTop 	,Config.Etapas);
	DataNxN.St.Bat.sTop	 = 	replicateMat3_3(	DataNxN.St.Bat.sTop 	,Config.Etapas);
	DataNxN.St.Bat.xiTop	 = 	replicateMat3_3(	DataNxN.St.Bat.xiTop 	,Config.Etapas);
	DataNxN.St.Bat.wOm	 = 	replicateMat3_3(	DataNxN.St.Bat.wOm 	,Config.Etapas);
	DataNxN.St.Bat.wU	 = 	replicateMat3_3(	DataNxN.St.Bat.wU 	,Config.Etapas);

	DataNxN.temp	 = 	replicateMat3_1(	DataNxN.temp(Et)	,Buses);

	% DataNxN.Util.aT	 = 	replicateMat4_3_4(	DataNxN.Util.aT	,Config.Etapas,2);
	% DataNxN.Util.aE	 = 	replicateMat3_3(	DataNxN.Util.aE	,Config.Etapas);
	DataNxN.Util.betaE	 = 	replicateMat3_2(	DataNxN.Util.betaE(:,Et)	,1);
	DataNxN.Util.betaT	 = 	replicateMat4_3_4(	DataNxN.Util.betaT(:,Et)	,2);
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

function [A] = replicateMat4_3_4(V,D4)
	A = zeros(size(V,1),1,size(V,2),D4);
    for d4 = 1:D4
        A(:,1,:,d4) = V;
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
