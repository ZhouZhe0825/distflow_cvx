function printPv(Header, Var, Data, outFilename)

	if isfield(Var, 'Gen')
		if isfield(Var.Gen, 'Pv')
			nodPv = find(matOverTime(Data.Gen.Pv.I) == 1);

			for i = 1:length(nodPv)
				nod = nodPv(i);
				pPv = Var.Gen.Pv.pPv(nod, 1, :);
				qPv = Var.Gen.Pv.qPv(nod, 1, :);
				s = Var.Gen.Pv.s(nod, 1, :);
				h_s = sqrt(pPv.^2 + qPv.^2)./s;
				cv_s = Data.Gen.Pv.cv(nod, 1, :).*s;
				xi = Var.Gen.Pv.xi(nod, 1, :);
				h_xi = (pPv.^2 + qPv.^2)./xi;
				cr_xi = Data.Gen.Pv.cr(nod, 1, :).*xi;
				Pv = [pPv;qPv;s;h_s;cv_s;xi;h_xi;cr_xi];



				
				rowHeader = cell(8,1);
				rowHeader{1} = 'pPv';
				rowHeader{2} = 'qPv';
				rowHeader{3} = 's';
				rowHeader{4} = 'Holgura s';
				rowHeader{5} = 'cv s';
				rowHeader{6} = 'x';
				rowHeader{7} = 'Holgura x';
				rowHeader{8} = 'cv s';
				sheetName = ['Pv_' num2str(nod)];
				printVarNx1xT(Pv, rowHeader, Header, outFilename, sheetName);
			end
		end
	end
end
