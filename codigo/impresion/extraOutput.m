function extraOutput(Var, Data, Config, Header, outFilename)

    n = size(Data.Red.Branch.T,1);
	TotalT = matOverTime(Data.Red.Branch.T);

    G = find(Data.Gen.Tras.I == 1);
    NnoG = setdiff((1:n), G)';
    TSalientesG = Data.Red.Branch.T;
    TEntrantesG = Data.Red.Branch.T;
    TnoG = Data.Red.Branch.T;
    TnoG(G,:,:) = zeros(length(G),n,Config.Etapas);
    TSalientesG = TSalientesG - TnoG;
    TnoG(:,G,:) = zeros(n,length(G),Config.Etapas);

    P = Var.Red.Branch.P;
    Q = Var.Red.Branch.Q;
    l = Var.Red.Branch.l;
    z = Var.Red.Branch.z;
    nv = Var.Red.Branch.nv;
    v = Var.Red.Bus.v;

    tvEqL = (nv - repmat(permute(v, [2 1 3]), [n 1 1])).*Data.Red.Branch.T;
    tvEqR = - 2 * (Data.Red.Branch.r .* P + Data.Red.Branch.x .* Q) + (Data.Red.Branch.r.^2 + Data.Red.Branch.x.^2) .* l;
    tvEqz = repmat(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2, [1 n 1]).*(1-z);
    l0 = (tvEqL + tvEqR - tvEqz).*TnoG;
    G0 = (tvEqL + tvEqR + tvEqz).*TnoG;
    
    h_p = abs((Var.Red.Bus.pC - Var.Red.Bus.pG) ./ Var.Red.Bus.pN);
    h_q = abs((Var.Red.Bus.qC - Var.Red.Bus.qG) ./ Var.Red.Bus.qN);
    
	lzr = l .* Data.Red.Branch.r .* z;

	PQv = (P.^2 + Q.^2) ./ repmat(v, [1 n 1]).*Data.Red.Branch.T;
	lRel = abs(PQv) ./ (abs(l)+eps);
	lRelz = lRel .* z;

    indl = intersect(intersect(find(lzr > 1e-8),find(lRelz < .9)),find(lRelz > 1e-5));
    lRelz_prob = lRelz*0;
    lRelz_prob(indl) = lRelz(indl);
   
    printVarNx1xT(h_p, Header.Bus, Header.Main, outFilename, 'pC - pG div pN');
    printVarNx1xT(h_q, Header.Bus, Header.Main, outFilename, 'qC - qG div qN');

	printVarNx1xT(TreeMatTimToVectTim(G0,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'EqV >= 0');
	printVarNx1xT(TreeMatTimToVectTim(l0,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'EqV <= 0');

    printVarNx1xT(TreeMatTimToVectTim(PQv,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'P2+Q2 div v');
	printVarNx1xT(TreeMatTimToVectTim(lzr,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'l.z.r');
	printVarNx1xT(TreeMatTimToVectTim(lRel,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'Holgura l');
	printVarNx1xT(TreeMatTimToVectTim(lRelz,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'Holgura l.z');
	printVarNx1xT(TreeMatTimToVectTim(lRelz_prob,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'Holguras l prob');
	
	if isfield(Var, 'Dual')
% 		dPn = zeros(n,Config.Etapas);
% 		dQn = zeros(n,Config.Etapas);
        dPn = squeeze(Var.Dual.dPn);
        dQn = squeeze(Var.Dual.dQn);
		printVarNx1xT(dPn, Header.Bus, Header.Main, outFilename, 'dPn');
		printVarNx1xT(dQn, Header.Bus, Header.Main, outFilename, 'dQn');
    end
	
	
end
