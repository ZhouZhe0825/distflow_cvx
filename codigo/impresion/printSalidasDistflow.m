function printSalidasDistflow(Var, Data, Config, cantTaps, cantCaps, cantCarg, outFilename, optEv, muEv, lambdaEv, DifPEv, DifQEv)

%     n = size(Data.Red.Branch.T,1);
%     nodos = (1:n)';

%     Output.opt(1) = sum(opt(1,:));


    TotalT = matOverTime(Data.Red.Branch.T);

%     [rowT, colT, ~] = find(TotalT == 1);


	Header = createHeader(Var, Data, Config, cantTaps, cantCaps, cantCarg);

	printVarNx1xT(sqrt(Var.Red.Bus.v), Header.Bus, Header.Main, outFilename, 'U');
	printVarNx1xT(TreeMatTimToVectTim(Var.Red.Branch.l,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'l');
	printVarNx1xT(TreeMatTimToVectTim(Var.Red.Branch.P,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'P');
	printVarNx1xT(TreeMatTimToVectTim(Var.Red.Branch.Q,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'Q');
	printVarNx1xT(TreeMatTimToVectTim(Var.Red.Branch.z,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'z');
	printVarNx1xT(TreeMatTimToVectTim(Var.Red.Branch.y,Config.Etapas,TotalT), Header.Branch, Header.Main, outFilename, 'y');
	printVarNx1xT(Var.Red.Bus.w, Header.Bus, Header.Main, outFilename, 'w');
	printVarNx1xT(Var.Red.Bus.cDv, Header.Bus, Header.Main, outFilename, 'cDv');
	printVarNx1xT(Var.Red.Bus.nn, Header.Bus, Header.Main, outFilename, 'nn');
	printVarNx1xT(Var.Red.Bus.nv, Header.Bus, Header.Main, outFilename, 'nv');
	printVarNx1xT(Var.Red.Bus.pG, Header.Bus, Header.Main, outFilename, 'pG');
	printVarNx1xT(Var.Red.Bus.pN, Header.Bus, Header.Main, outFilename, 'pN');
	printVarNx1xT(Var.Red.Bus.pC, Header.Bus, Header.Main, outFilename, 'pC');
	printVarNx1xT(Var.Red.Bus.qG, Header.Bus, Header.Main, outFilename, 'qG');
	printVarNx1xT(Var.Red.Bus.qN, Header.Bus, Header.Main, outFilename, 'qN');
	printVarNx1xT(Var.Red.Bus.qC, Header.Bus, Header.Main, outFilename, 'qC');
	
	printVarNx1xT(Var.ClRes.pC, Header.Bus, Header.Main, outFilename, 'pClRes');
	printVarNx1xT(Var.ClRes.qC, Header.Bus, Header.Main, outFilename, 'qClRes');
	if isfield(Var, 'ClRes')
		if isfield(Var.ClRes, 'Tvar')
			printVarNx1xT(Var.ClRes.Tvar, Header.Bus, Header.Main, outFilename, 'Tvar');
		end
	end
	printVarNx1xT(squeeze(Var.ClRes.pCApp(:,:,:,1)), Header.Bus, Header.Main, outFilename, 'pCApp_1');
	printVarNx1xT(squeeze(Var.ClRes.pCApp(:,:,:,2)), Header.Bus, Header.Main, outFilename, 'pCApp_2');
	printVarNx1xT(Var.Red.Bus.pN, Header.Bus, Header.Main, outFilename, 'pNRed');
	printVarNx1xT(Var.Red.Bus.qN, Header.Bus, Header.Main, outFilename, 'qNRed');
	printVarNx1xT(Var.Red.Bus.qCp, Header.Bus, Header.Main, outFilename, 'qCp');
	printVarNx1xT(Var.Red.Bus.PTras, Header.Bus, Header.Main, outFilename, 'Ptras');
	printVarNx1xT(Var.Red.Bus.QTras, Header.Bus, Header.Main, outFilename, 'Qtras');

	printClNI(Header.Main, Var, Data, outFilename);
	printPv(Header.Main, Var, Data, outFilename);
	printBat(Header.Main, Var, Data, outFilename);
	printDfig(Header.Main, Var, Data, outFilename);



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
    extraOutput(Var, Data, Config, Header, outFilename);
	
end
