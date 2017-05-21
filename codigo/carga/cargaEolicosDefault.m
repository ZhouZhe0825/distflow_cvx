function [Data] = cargaEolicosDefault(Data,Eolicos)

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
    
    for i = 1:length(Eolicos)
        [vars] = loadCsvData(Eolicos(i).fileG,n);

        P_mec = [];
        n_ = [];
        
        lenVar = length(vars);
        if lenVar == 2
            for j = 1:lenVar
                if strcmp(vars(j).name, 'P_mec') && vars(j).undefBus
                    P_mec = vars(j).data;
                elseif  strcmp(vars(j).name, 'n_') && vars(j).undefBus
                    n_ = vars(j).data;
                end
            end
        end
        if ~isempty(P_mec) && ~isempty(n_)
            [w,ws,Np,G,R_,lambda_opt,Omega,rho,N_er,RW_12,LW_12,RW_13,LW_13,RW_45,LW_45,...
                lW_12_Top,lW_13_Top,lW_45_Top,sW_3_Top,sW_5_Top,vmm,vm,vM,vMM,cvW_3,cvW_5,crW_3,crW_5,...
                P_nMec,c1,c2,c3,c4,c5,c6,c7,uLow_2,uTop_2,uLow_3,uTop_3,PQnorm_2,PQnorm_3,C_plb,vv] = Eolicos(i).type();

            Data.Gen.DFIG.P_mec(Eolicos(i).nod,:) = P_mec;
            Data.Gen.DFIG.n_(Eolicos(i).nod,:) = n_;
 
            Data.Gen.DFIG.N_er(Eolicos(i).nod) = N_er;
            Data.Gen.DFIG.w(Eolicos(i).nod) = w;
            Data.Gen.DFIG.rho(Eolicos(i).nod) = rho;
            Data.Gen.DFIG.R_(Eolicos(i).nod) = R_;
            Data.Gen.DFIG.C_plb(Eolicos(i).nod) = C_plb;
            Data.Gen.DFIG.vv(Eolicos(i).nod) = vv;

            Data.Gen.DFIG.vmm(Eolicos(i).nod) = vmm;
            Data.Gen.DFIG.vm(Eolicos(i).nod) = vm;
            Data.Gen.DFIG.vM(Eolicos(i).nod) = vM;
            Data.Gen.DFIG.vMM(Eolicos(i).nod) = vMM;

            Data.Gen.DFIG.P_nMec(Eolicos(i).nod) = P_nMec;
            Data.Gen.DFIG.c1(Eolicos(i).nod) = c1;
            Data.Gen.DFIG.c2(Eolicos(i).nod) = c2;
            Data.Gen.DFIG.c3(Eolicos(i).nod) = c3;
            Data.Gen.DFIG.c4(Eolicos(i).nod) = c4;
            Data.Gen.DFIG.c5(Eolicos(i).nod) = c5;
            Data.Gen.DFIG.c6(Eolicos(i).nod) = c6;
            Data.Gen.DFIG.c7(Eolicos(i).nod) = c7;

            Data.Gen.DFIG.Omega(Eolicos(i).nod) = Omega;
            Data.Gen.DFIG.G(Eolicos(i).nod) = G;
            Data.Gen.DFIG.Np(Eolicos(i).nod) = Np;
            Data.Gen.DFIG.rho(Eolicos(i).nod) = rho;
            Data.Gen.DFIG.ws(Eolicos(i).nod) = ws;
            Data.Gen.DFIG.lambda_opt(Eolicos(i).nod) = lambda_opt;

            Data.Gen.DFIG.rIE(Eolicos(i).nod,:) = RW_12;
            Data.Gen.DFIG.rIF(Eolicos(i).nod,:) = RW_13;
            Data.Gen.DFIG.rOR(Eolicos(i).nod,:) = RW_45;

            Data.Gen.DFIG.xIE(Eolicos(i).nod,:) = repmat(Data.Gen.DFIG.w(Eolicos(i).nod)*LW_12, [1,et]);
            Data.Gen.DFIG.xIF(Eolicos(i).nod,:) = repmat(Data.Gen.DFIG.w(Eolicos(i).nod)*LW_13, [1,et]);
            Data.Gen.DFIG.xOR(Eolicos(i).nod,:) = repmat(Data.Gen.DFIG.w(Eolicos(i).nod)*LW_45, [1,et]);

            Data.Gen.DFIG.lTopIE(Eolicos(i).nod,:) = lW_12_Top;
            Data.Gen.DFIG.lTopIF(Eolicos(i).nod,:) = lW_13_Top;
            Data.Gen.DFIG.lTopOR(Eolicos(i).nod,:) = lW_45_Top;

            Data.Gen.DFIG.cvF(Eolicos(i).nod,:) = cvW_3;
            Data.Gen.DFIG.cvR(Eolicos(i).nod,:) = cvW_5;

            Data.Gen.DFIG.crF(Eolicos(i).nod,:) = crW_3;
            Data.Gen.DFIG.crR(Eolicos(i).nod,:) = crW_5;

            Data.Gen.DFIG.uLowE(Eolicos(i).nod,:) = uLow_2;
            Data.Gen.DFIG.uLowF(Eolicos(i).nod,:) = uLow_3;

            Data.Gen.DFIG.uTopE(Eolicos(i).nod,:) = uTop_2;
            Data.Gen.DFIG.uTopF(Eolicos(i).nod,:) = uTop_3;

            Data.Gen.DFIG.sTopF(Eolicos(i).nod,:) = sW_3_Top;
            Data.Gen.DFIG.sTopR(Eolicos(i).nod,:) = sW_5_Top;

            Data.Gen.DFIG.xiTopF(Eolicos(i).nod,:) = sW_3_Top^2;
            Data.Gen.DFIG.xiTopR(Eolicos(i).nod,:) = sW_5_Top^2;

            Data.Gen.DFIG.I(Eolicos(i).nod,:) = 1;

            Data.Gen.DFIG.PQnormIE(Eolicos(i).nod,:) = PQnorm_2;
            Data.Gen.DFIG.PQnormIF(Eolicos(i).nod,:) = PQnorm_3;
        end
    end


    