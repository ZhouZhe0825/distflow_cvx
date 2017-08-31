function extraOutput(Var, Data, Config, Header, outFilename)

    n = size(Data.Red.Branch.T,1);

	VertI = VertIMat(Data.Red.Branch.T);
	VertJ = VertJMat(Data.Red.Branch.T);
	
	
    G = find(Data.Gen.Tras.I == 1);
%     NnoG = setdiff((1:n), G)';
%     TSalientesG = Data.Red.Branch.T;
%     TEntrantesG = Data.Red.Branch.T;
%     TnoG = Data.Red.Branch.T;
%     TnoG(G,:,:) = zeros(length(G),n,Config.Etapas);
%     TSalientesG = TSalientesG - TnoG;
%     TnoG(:,G,:) = zeros(n,length(G),Config.Etapas);

    P = Var.Red.Branch.P;
    Q = Var.Red.Branch.Q;
    l = Var.Red.Branch.l;
    z = Var.Red.Branch.z;
    nv = Var.Red.Branch.nv;
    v = Var.Red.Bus.v;

    tvEqL = nv - VertJ*v;
    tvEqR = - 2 * (Data.Red.Branch.r .* P + Data.Red.Branch.x .* Q) + (Data.Red.Branch.r.^2 + Data.Red.Branch.x.^2) .* l;
    tvEqz = (VertI*(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2)).*(1-z);
%     l0 = (tvEqL + tvEqR - tvEqz).*TnoG;
%     G0 = (tvEqL + tvEqR + tvEqz).*TnoG;
    l0 = (tvEqL + tvEqR - tvEqz);
    G0 = (tvEqL + tvEqR + tvEqz);
    
    h_p = abs((Var.Red.Bus.pC - Var.Red.Bus.pG) ./ Var.Red.Bus.pN);
    h_q = abs((Var.Red.Bus.qC - Var.Red.Bus.qG) ./ Var.Red.Bus.qN);
    
	lzr = l .* Data.Red.Branch.r .* z;

	PQv = (P.^2 + Q.^2) ./ (VertI*v);
	lRel = abs(PQv) ./ (abs(l)+eps);
	lRelz = lRel .* z;

    indl = intersect(intersect(find(lzr > 1e-8),find(lRelz < .9)),find(lRelz > 1e-5));
    lRelz_prob = lRelz*0;
    lRelz_prob(indl) = lRelz(indl);
   
    printVarMxT(h_p, Header.Bus, Header.Main, outFilename, 'pC - pG div pN');
    printVarMxT(h_q, Header.Bus, Header.Main, outFilename, 'qC - qG div qN');

	printVarMxT(G0, Header.Branch, Header.Main, outFilename, 'EqV >= 0');
	printVarMxT(l0, Header.Branch, Header.Main, outFilename, 'EqV <= 0');

    printVarMxT(PQv, Header.Branch, Header.Main, outFilename, 'P2+Q2 div v');
	printVarMxT(lzr, Header.Branch, Header.Main, outFilename, 'l.z.r');
	printVarMxT(lRel, Header.Branch, Header.Main, outFilename, 'Holgura l');
	printVarMxT(lRelz, Header.Branch, Header.Main, outFilename, 'Holgura l.z');
	printVarMxT(lRelz_prob, Header.Branch, Header.Main, outFilename, 'Holguras l prob');
	
	if isfield(Var, 'Dual')
        dPn = Var.Dual.dPn;
        dQn = Var.Dual.dQn;
		printVarMxT(dPn, Header.Bus, Header.Main, outFilename, 'dPn');
		printVarMxT(dQn, Header.Bus, Header.Main, outFilename, 'dQn');
    end
	
	
end
