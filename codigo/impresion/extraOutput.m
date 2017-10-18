function extraOutput(Var, Data, Config, Header, outFilename)

%     n = size(Data.Red.Branch.T,1);

	VertI = VertIMat(Data.Red.Branch.T);
	VertJ = VertJMat(Data.Red.Branch.T);
	
%     G = find(Data.Gen.Tras.I == 1);

    P = Var.Red.Branch.P;
    Q = Var.Red.Branch.Q;
    l = Var.Red.Branch.l;
    z = Var.Red.Branch.z;
    nv = Var.Red.Branch.nv;
    v = Var.Red.Bus.v;
    pG = Var.Red.Bus.pG;
    qG = Var.Red.Bus.qG;
    pC = Var.Red.Bus.pC;
    qC = Var.Red.Bus.qC;
    pN = Var.Red.Bus.pN;
    qN = Var.Red.Bus.qN;

    tvEqL = nv - VertJ*v;
    tvEqR = - 2 * (Data.Red.Branch.r .* P + Data.Red.Branch.x .* Q) + (Data.Red.Branch.r.^2 + Data.Red.Branch.x.^2) .* l;
    tvEqz = (VertI*(Data.Red.Bus.uTop.^2 - Data.Red.Bus.uLow.^2)).*(1-z);
    l0 = (tvEqL + tvEqR - tvEqz);
    G0 = (tvEqL + tvEqR + tvEqz);
    
    h_p_rel = abs((pC - pG) ./ pN);
    h_q_rel = abs((qC - qG) ./ qN);
    h_p_abs = abs(pC - pG - pN);
    h_q_abs = abs(qC - qG - qN);
    
	lzr = l .* Data.Red.Branch.r .* z;

	PQv = (P.^2 + Q.^2) ./ (VertI*v);
	lRel = abs(PQv) ./ (abs(l)+eps);
	lRelz = lRel .* z;

    indl = intersect(intersect(find(lzr > 1e-8),find(lRelz < .9)),find(lRelz > 1e-5));
    lRelz_prob = lRelz*0;
    lRelz_prob(indl) = lRelz(indl);
   
    printVarMxT(h_p_rel, Header.Bus, Header.Main, outFilename, 'pC - pG div pN');
    printVarMxT(h_q_rel, Header.Bus, Header.Main, outFilename, 'qC - qG div qN');

    printVarMxT(h_p_abs, Header.Bus, Header.Main, outFilename, 'pC - pG - pN');
    printVarMxT(h_q_abs, Header.Bus, Header.Main, outFilename, 'qC - qG - qN');

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
