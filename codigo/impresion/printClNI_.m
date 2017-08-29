function printClNI_(Header, Var, Data, outFilename)

	if isfield(Var, 'ClNI')
		nodCh = find(Data.ClNI.I == 1);

		for i = 1:length(nodCh)
			nod = nodCh(i);
			pC = Var.ClNI.pC(nod, :);
			qC = Var.ClNI.qC(nod, :);
			on = Var.ClNI.on(nod, :);
			start = Var.ClNI.start(nod, :);
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
