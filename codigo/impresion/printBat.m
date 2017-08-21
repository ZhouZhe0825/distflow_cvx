function printBat(Header, Var, Data, outFilename)

	if isfield(Var, 'St')
		if isfield(Var.St, 'Bat')
			nodBat = find(matOverTime(Data.St.Bat.I) == 1);

			for i = 1:length(nodBat)
				nod = nodBat(i);
				pStb = Var.St.Bat.pStb(nod, 1, :);
				pStgbC = Var.St.Bat.pStgbC(nod, 1, :);
				pStgbD = Var.St.Bat.pStgbD(nod, 1, :);
				qStb = Var.St.Bat.qStb(nod, 1, :);
				EStb = Var.St.Bat.EStb(nod, 1, :);
				pStgb = Var.St.Bat.pStgb(nod, 1, :);
				sStb = Var.St.Bat.sStb(nod, 1, :);
				h_sStb = sqrt(pStgb.^2 + qStb.^2)./sStb;
				cv_sStb = Data.St.Bat.cv(nod,1,:).*sStb;
				xiStb = Var.St.Bat.xiStb(nod, 1, :);
				h_xiStb = (pStgb.^2 + qStb.^2)./xiStb;
				cr_xiStb = Data.St.Bat.cr(nod,1,:).*xiStb;

				rowHeader = cell(12,1);
				rowHeader{1} = 'pStb';
				rowHeader{2} = 'pStgbC';
				rowHeader{3} = 'pStgbD';
				rowHeader{4} = 'qStb';
				rowHeader{5} = 'EStb';
				rowHeader{6} = 'pStgb';
				rowHeader{7} = 'sStb';
				rowHeader{8} = 'Holgura sStb';
				rowHeader{9} = 'cv sStb';
				rowHeader{10} = 'xiStb';
				rowHeader{11} = 'Holgura xiStb';
				rowHeader{12} = 'cr xiStb';

				Bat = [pStb; pStgbC; pStgbD; qStb; EStb; pStgb; sStb; h_sStb; cv_sStb; xiStb; h_xiStb; cr_xiStb];
				sheetName = ['Bat_' num2str(nod)];
				printVarNx1xT(Bat, rowHeader, Header, outFilename, sheetName);
			end
		end
	end
end
