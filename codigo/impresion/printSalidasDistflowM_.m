function printSalidasDistflowM_(Var, Data, Config, outFilename, optEv, muEv, lambdaEv, DifPEv, DifQEv)

%     n = size(Data.Red.Branch.T,1);
%     nodos = (1:n)';

%     Output.opt(1) = sum(opt(1,:));



	Header = createHeader_(Var, Data, Config);

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
	
	printVarMxT(Var.ClRes.pC, Header.Bus, Header.Main, outFilename, 'pClRes');
	printVarMxT(Var.ClRes.qC, Header.Bus, Header.Main, outFilename, 'qClRes');
	if isfield(Var, 'ClRes')
		if isfield(Var.ClRes, 'Tvar')
			printVarMxT(Var.ClRes.Tvar, Header.Bus, Header.Main, outFilename, 'Tvar');
		end
	end
	printVarMxT(squeeze(Var.ClRes.pCApp(:,:,1)), Header.Bus, Header.Main, outFilename, 'pCApp_1');
	printVarMxT(squeeze(Var.ClRes.pCApp(:,:,2)), Header.Bus, Header.Main, outFilename, 'pCApp_2');
	printVarMxT(squeeze(Var.ClRes.qCApp(:,:,1)), Header.Bus, Header.Main, outFilename, 'qCApp_1');
	printVarMxT(squeeze(Var.ClRes.qCApp(:,:,2)), Header.Bus, Header.Main, outFilename, 'qCApp_2');
	printVarMxT(Var.Red.Bus.pN, Header.Bus, Header.Main, outFilename, 'pNRed');
	printVarMxT(Var.Red.Bus.qN, Header.Bus, Header.Main, outFilename, 'qNRed');
	printVarMxT(Var.Red.Bus.qCp, Header.Bus, Header.Main, outFilename, 'qCp');
	printVarMxT(Var.Red.Bus.PTras, Header.Bus, Header.Main, outFilename, 'Ptras');
	printVarMxT(Var.Red.Bus.QTras, Header.Bus, Header.Main, outFilename, 'Qtras');

    printCost_(Header.Main, Var, Data, outFilename);
    printClNI_(Header.Main, Var, Data, outFilename);
	printBasic_(Header.Main, Var, Data, outFilename);
    printPv_(Header.Main, Var, Data, outFilename);
	printBat_(Header.Main, Var, Data, outFilename);
	printDfig_(Header.Main, Var, Data, outFilename);



    % if isfield(Data, 'Fixed')
        % dPnPrint = [Header.Main;num2cell([nodos Output.dPn])];
        % dQnPrint = [Header.Main;num2cell([nodos Output.dQn])];

        % xlswrite([outFilename, '.xlsx'], dPnPrint, 'dPn');
        % xlswrite([outFilename, '.xlsx'], dQnPrint, 'dQn');
    %     end


    it = min([size(optEv,1), size(muEv, 3), size(lambdaEv, 3), size(DifPEv, 3), size(DifQEv, 3)]);

    if it > 0
        xlswrite([outFilename, '.xlsx'],reshape(optEv, [size(optEv,1) size(optEv,2)]), 'optEv');
        muEvOut = zeros(size(muEv,1)*size(muEv,3)+it-1, size(muEv,2));
        lambdaEvOut = zeros(size(lambdaEv,1)*size(lambdaEv,3)+it-1, size(lambdaEv,2));
        DifPEvOut = zeros(size(DifPEv,1)*size(DifPEv,3)+it-1, size(DifPEv,2));
        DifQEvOut = zeros(size(DifQEv,1)*size(DifQEv,3)+it-1, size(DifQEv,2));
        for i = 1:it
            lambdaEvOut((size(lambdaEv,1)*(i-1)+i:size(lambdaEv,1)*i+i-1),:) = reshape(lambdaEv(:,:,i), [size(lambdaEv,1) size(lambdaEv,2)]);
            DifPEvOut((size(DifPEv,1)*(i-1)+i:size(DifPEv,1)*i+i-1),:) = reshape(DifPEv(:,:,i), [size(DifPEv,1) size(DifPEv,2)]);
            DifQEvOut((size(DifQEv,1)*(i-1)+i:size(DifQEv,1)*i+i-1),:) = reshape(DifQEv(:,:,i), [size(DifQEv,1) size(DifQEv,2)]);
            muEvOut((size(muEv,1)*(i-1)+i:size(muEv,1)*i+i-1),:) = reshape(muEv(:,:,i), [size(muEv,1) size(muEv,2)]);
        end
        xlswrite([outFilename, '.xlsx'],muEvOut, 'muEv');
        xlswrite([outFilename, '.xlsx'],lambdaEvOut, 'lambdaEv');
        xlswrite([outFilename, '.xlsx'],DifPEvOut, 'DifPEv');
        xlswrite([outFilename, '.xlsx'],DifQEvOut, 'DifQEv');
    end
    extraOutput_(Var, Data, Config, Header, outFilename);
	
end
