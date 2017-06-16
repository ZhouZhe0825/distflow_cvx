function printSalidasDistflowM(Var, Data, Config, outFilename, optEv, muEv, lambdaEv, DifPEv, DifQEv)

    TotalT = matOverTime(Data.Red.Branch.T);
    
    D.Red.Branch.T = TotalT;

    [Var] = VarM2NxN(Var,D);

	printSalidasDistflowNxN(Var, Data, Config, outFilename, optEv, muEv, lambdaEv, DifPEv, DifQEv);
end
