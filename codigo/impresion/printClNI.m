function printClNI(Header, Var, Data, outFilename)

	if isfield(Var, 'ClNI')
		nodCh = find(matOverTime(Data.ClNI.I) == 1);

		for i = 1:length(nodCh)
			nod = nodCh(i);
			pC = Var.ClNI.pC(nod, 1, :);
			qC = Var.ClNI.qC(nod, 1, :);
			on = Var.ClNI.on(nod, 1, :);
			start = Var.ClNI.start(nod, 1, :);
			ClNI = [pC;qC;on;start];
			rowHeader = cell(4,1);
			rowHeader{1} = 'pC';
			rowHeader{2} = 'qC';
			rowHeader{3} = 'on';
			rowHeader{4} = 'start';
			sheetName = ['ClNI_' num2str(nod)];
			printVarNx1xT(ClNI, rowHeader, Header, outFilename, sheetName);
		end
	end
end
