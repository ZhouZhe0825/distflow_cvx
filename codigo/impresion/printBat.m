function printBat(Header, Var, Data, outFilename)

	if isfield(Var, 'St')
		if isfield(Var.St, 'Bat')
			nodBat = find(matOverTime(Data.St.Bat.I) == 1);

			for i = 1:length(nodBat)
				nod = nodBat(i);
				pStb = Var.St.Bat.pStb(nod, 1, :);
				qStb = Var.St.Bat.qStb(nod, 1, :);
				EStb = Var.St.Bat.EStb(nod, 1, :);
				pStgb = Var.St.Bat.pStgb(nod, 1, :);
				sStb = Var.St.Bat.sStb(nod, 1, :);
				xiStb = Var.St.Bat.xiStb(nod, 1, :);

				rowHeader = cell(6,1);
				rowHeader{1} = 'pStb';
				rowHeader{2} = 'qStb';
				rowHeader{3} = 'EStb';
				rowHeader{4} = 'pStgb';
				rowHeader{5} = 'sStb';
				rowHeader{6} = 'xiStb';

				Bat = [pStb; qStb; EStb; pStgb; sStb; xiStb];
				sheetName = ['Bat_' num2str(nod)];
				printVarNx1xT(Bat, rowHeader, Header, outFilename, sheetName);
			end
		end
	end
end
