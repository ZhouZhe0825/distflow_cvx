function printPv(Header, Var, Data, outFilename)

	if isfield(Var, 'Gen')
		if isfield(Var.Gen, 'Pv')
			nodPv = find(matOverTime(Data.Gen.Pv.I) == 1);

			for i = 1:length(nodPv)
				nod = nodPv(i);
				pPv = Var.Gen.Pv.pPv(nod, 1, :);
				qPv = Var.Gen.Pv.qPv(nod, 1, :);
				sPv = Var.Gen.Pv.s(nod, 1, :);
				h_sPv = sqrt(pPv.^2 + qPv.^2)./sPv;
				cv_sPv = Data.Gen.Pv.cv(nod, 1, :).*sPv;
				xiPv = Var.Gen.Pv.xi(nod, 1, :);
				h_xiPv = (pPv.^2 + qPv.^2)./xiPv;
				cr_xiPv = Data.Gen.Pv.cr(nod, 1, :).*xiPv;
				Pv = [pPv;qPv;sPv;h_sPv;cv_sPv;xiPv;h_xiPv;cr_xiPv];



				
				rowHeader = cell(8,1);
				rowHeader{1} = 'pPv';
				rowHeader{2} = 'qPv';
				rowHeader{3} = 'sPv';
				rowHeader{4} = 'Holgura sPv';
				rowHeader{5} = 'cv sPv';
				rowHeader{6} = 'xiPv';
				rowHeader{7} = 'Holgura xiPv';
				rowHeader{8} = 'cr xiPv';
				sheetName = ['Pv_' num2str(nod)];
				printVarNx1xT(Pv, rowHeader, Header, outFilename, sheetName);
			end
		end
	end
end
