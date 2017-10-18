function printClRes(Header, Var, Data, outFilename)

	if isfield(Var, 'ClRes')
		if isfield(Var.Gen, 'Basic')
            printVarMxT(Var.ClRes.pC, Header.Bus, Header.Main, outFilename, 'pClRes');
            printVarMxT(Var.ClRes.qC, Header.Bus, Header.Main, outFilename, 'qClRes');
            if isfield(Var.ClRes, 'Tvar')
                printVarMxT(Var.ClRes.Tvar, Header.Bus, Header.Main, outFilename, 'Tvar');
            end
            for i=1:size(Var.ClRes.pCApp,3)
                printVarMxT(squeeze(Var.ClRes.pCApp(:,:,i)), Header.Bus, Header.Main, outFilename, ['pCApp_', num2str(i)]);
                printVarMxT(squeeze(Var.ClRes.qCApp(:,:,i)), Header.Bus, Header.Main, outFilename, ['qCApp_', num2str(i)]);
            end
		end
	end
end
