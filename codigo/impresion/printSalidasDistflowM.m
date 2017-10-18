function printSalidasDistflowM(Var, Data, Config, outFilename, optEv, muEv, lambdaEv, DifPEv, DifQEv)

%     n = size(Data.Red.Branch.T,1);
%     nodos = (1:n)';

%     Output.opt(1) = sum(opt(1,:));



	Header = createHeader(Var, Data, Config);

	printVarMxT(sqrt(Var.Red.Bus.v), Header.Bus, Header.Main, outFilename, 'U');
	printVarMxT(Var.Red.Branch.l, Header.Branch, Header.Main, outFilename, 'l');
	printVarMxT(Var.Red.Branch.P, Header.Branch, Header.Main, outFilename, 'P');
	printVarMxT(Var.Red.Branch.Q, Header.Branch, Header.Main, outFilename, 'Q');
	printVarMxT(Var.Red.Branch.z, Header.Branch, Header.Main, outFilename, 'z');
	printVarMxT(Var.Red.Branch.y, Header.Branch, Header.Main, outFilename, 'y');
	printVarMxT(Var.Red.Bus.w, Header.Bus, Header.Main, outFilename, 'w');
	printVarMxT(Var.Red.Bus.cDv, Header.Bus, Header.Main, outFilename, 'cDv');
	printVarMxT(Var.Red.Branch.nn, Header.Branch, Header.Main, outFilename, 'nn');
	printVarMxT(Var.Red.Branch.nv, Header.Branch, Header.Main, outFilename, 'nv');
	printVarMxT(Var.Red.Bus.pG, Header.Bus, Header.Main, outFilename, 'pG');
	printVarMxT(Var.Red.Bus.pN, Header.Bus, Header.Main, outFilename, 'pN');
	printVarMxT(Var.Red.Bus.pC, Header.Bus, Header.Main, outFilename, 'pC');
	printVarMxT(Var.Red.Bus.qG, Header.Bus, Header.Main, outFilename, 'qG');
	printVarMxT(Var.Red.Bus.qN, Header.Bus, Header.Main, outFilename, 'qN');
	printVarMxT(Var.Red.Bus.qC, Header.Bus, Header.Main, outFilename, 'qC');
	printVarMxT(Var.Red.Bus.qCp, Header.Bus, Header.Main, outFilename, 'qCp');

    printClRes(Header, Var, Data, outFilename);
    printCost(Header, Var, Data, outFilename);
    printClNI(Header, Var, Data, outFilename);
	printBasic(Header, Var, Data, outFilename);
    printPv(Header, Var, Data, outFilename);
	printBat(Header, Var, Data, outFilename);
	printDfig(Header, Var, Data, outFilename);

    extraOutput(Var, Data, Config, Header, outFilename);
	
end
