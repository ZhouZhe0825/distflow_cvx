function printBasic(Header, Var, Data, outFilename)

	if isfield(Var, 'Gen')
		if isfield(Var.Gen, 'Basic')
			nodGB = find(Data.Gen.Basic.I == 1);

			for i = 1:length(nodGB)
				nod = nodGB(i);
				pGB = Var.Gen.Basic.pGBas(nod, :);
				qGB = Var.Gen.Basic.qGBas(nod, :);
				GB = [pGB;qGB];

				
				rowHeader = cell(2,1);
				rowHeader{1} = 'p';
				rowHeader{2} = 'q';
				sheetName = ['GBas_' num2str(nod)];
				printVarNx1xT(GB, rowHeader, Header, outFilename, sheetName);
			end
		end
	end
end
