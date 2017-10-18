function printTras(Header, Var, Data, outFilename)

	if isfield(Var, 'Gen')
		if isfield(Var.Gen, 'Basic')
            printVarMxT(Var.Gen.Tras.pTras, Header.Bus, Header.Main, outFilename, 'pTras');
            printVarMxT(Var.Gen.Tras.qTras, Header.Bus, Header.Main, outFilename, 'qTras');
			end
		end
	end
end
