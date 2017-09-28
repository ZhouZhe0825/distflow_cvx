function [Data] = cargaEolicosDefault(Data,Eolicos)

    n = size(Data.Red.Branch.T,1);
    et = size(Data.Red.Bus.pCLow,2);

    for i = 1:length(Eolicos)
        [vars] = loadCsvDataSeries(Eolicos(i).fileG,n);

        P_mec = [];
        n_ = [];
        
        lenVar = length(vars);
        if lenVar == 2
            for j = 1:lenVar
                if strcmp(vars(j).name, 'P_mec') && vars(j).undefBus
                    P_mec = vars(j).data;
                elseif strcmp(vars(j).name, 'n_') && vars(j).undefBus
                    n_ = vars(j).data;
                end
            end
        end
        if ~isempty(P_mec) && ~isempty(n_)
            
            fields = { ...
                'w','ws','Np','G','R_','lambda_opt','Omega','rho','N_er', ...
                'RW_12','LW_12','RW_13','LW_13','RW_45','LW_45',...
                'lW_12_Top','lW_13_Top','lW_45_Top','sW_3_Top','sW_5_Top',...
                'vmm','vm','vM','vMM','cvW_3','cvW_5','crW_3','crW_5',...
                'P_nMec','c1','c2','c3','c4','c5','c6','c7','uLow_2',...
                'uTop_2','uLow_3','uTop_3','PQnorm_2','PQnorm_3','C_plb','vv',...
                'qWiLow','qWiTop'...
                };
            S = loadCsvDef(Eolicos(i).type(), fields);

            Data.Gen.DFIG.P_mec(Eolicos(i).nod,:) = P_mec;
            Data.Gen.DFIG.n_(Eolicos(i).nod,:) = n_;
 
            Data.Gen.DFIG.N_er(Eolicos(i).nod) = S.N_er;
            Data.Gen.DFIG.w(Eolicos(i).nod) = S.w;
            Data.Gen.DFIG.rho(Eolicos(i).nod) = S.rho;
            Data.Gen.DFIG.R_(Eolicos(i).nod) = S.R_;
            Data.Gen.DFIG.C_plb(Eolicos(i).nod) = S.C_plb;
            Data.Gen.DFIG.vv(Eolicos(i).nod) = S.vv;

            Data.Gen.DFIG.vmm(Eolicos(i).nod) = S.vmm;
            Data.Gen.DFIG.vm(Eolicos(i).nod) = S.vm;
            Data.Gen.DFIG.vM(Eolicos(i).nod) = S.vM;
            Data.Gen.DFIG.vMM(Eolicos(i).nod) = S.vMM;

            Data.Gen.DFIG.P_nMec(Eolicos(i).nod) = S.P_nMec;
            Data.Gen.DFIG.c1(Eolicos(i).nod) = S.c1;
            Data.Gen.DFIG.c2(Eolicos(i).nod) = S.c2;
            Data.Gen.DFIG.c3(Eolicos(i).nod) = S.c3;
            Data.Gen.DFIG.c4(Eolicos(i).nod) = S.c4;
            Data.Gen.DFIG.c5(Eolicos(i).nod) = S.c5;
            Data.Gen.DFIG.c6(Eolicos(i).nod) = S.c6;
            Data.Gen.DFIG.c7(Eolicos(i).nod) = S.c7;

            Data.Gen.DFIG.Omega(Eolicos(i).nod) = S.Omega;
            Data.Gen.DFIG.G(Eolicos(i).nod) = S.G;
            Data.Gen.DFIG.Np(Eolicos(i).nod) = S.Np;
            Data.Gen.DFIG.rho(Eolicos(i).nod) = S.rho;
            Data.Gen.DFIG.ws(Eolicos(i).nod) = S.ws;
            Data.Gen.DFIG.lambda_opt(Eolicos(i).nod) = S.lambda_opt;

            Data.Gen.DFIG.rIE(Eolicos(i).nod,:) = S.RW_12;
            Data.Gen.DFIG.rIF(Eolicos(i).nod,:) = S.RW_13;
            Data.Gen.DFIG.rOR(Eolicos(i).nod,:) = S.RW_45;

            Data.Gen.DFIG.xIE(Eolicos(i).nod,:) = repmat(Data.Gen.DFIG.w(Eolicos(i).nod)*S.LW_12, [1,et]);
            Data.Gen.DFIG.xIF(Eolicos(i).nod,:) = repmat(Data.Gen.DFIG.w(Eolicos(i).nod)*S.LW_13, [1,et]);
            Data.Gen.DFIG.xOR(Eolicos(i).nod,:) = repmat(Data.Gen.DFIG.w(Eolicos(i).nod)*S.LW_45, [1,et]);

            Data.Gen.DFIG.lTopIE(Eolicos(i).nod,:) = S.lW_12_Top;
            Data.Gen.DFIG.lTopIF(Eolicos(i).nod,:) = S.lW_13_Top;
            Data.Gen.DFIG.lTopOR(Eolicos(i).nod,:) = S.lW_45_Top;

            Data.Gen.DFIG.cvF(Eolicos(i).nod,:) = S.cvW_3;
            Data.Gen.DFIG.cvR(Eolicos(i).nod,:) = S.cvW_5;

            Data.Gen.DFIG.crF(Eolicos(i).nod,:) = S.crW_3;
            Data.Gen.DFIG.crR(Eolicos(i).nod,:) = S.crW_5;

            Data.Gen.DFIG.uLowE(Eolicos(i).nod,:) = S.uLow_2;
            Data.Gen.DFIG.uLowF(Eolicos(i).nod,:) = S.uLow_3;

            Data.Gen.DFIG.uTopE(Eolicos(i).nod,:) = S.uTop_2;
            Data.Gen.DFIG.uTopF(Eolicos(i).nod,:) = S.uTop_3;

            Data.Gen.DFIG.sTopF(Eolicos(i).nod,:) = S.sW_3_Top;
            Data.Gen.DFIG.sTopR(Eolicos(i).nod,:) = S.sW_5_Top;

            Data.Gen.DFIG.xiTopF(Eolicos(i).nod,:) = S.sW_3_Top^2;
            Data.Gen.DFIG.xiTopR(Eolicos(i).nod,:) = S.sW_5_Top^2;

            Data.Gen.DFIG.qWiLow(Eolicos(i).nod,:) = S.qWiLow;
            Data.Gen.DFIG.qWiTop(Eolicos(i).nod,:) = S.qWiTop;

            Data.Gen.DFIG.I(Eolicos(i).nod,:) = 1;

            Data.Gen.DFIG.PQnormIE(Eolicos(i).nod,:) = S.PQnorm_2;
            Data.Gen.DFIG.PQnormIF(Eolicos(i).nod,:) = S.PQnorm_3;

            Data.Gen.DFIG.pgIni(Eolicos(i).nod,:) = Eolicos(i).pgIni;
            Data.Gen.DFIG.qgIni(Eolicos(i).nod,:) = Eolicos(i).qgIni;
        end
    end


    