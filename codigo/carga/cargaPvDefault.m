function [Data] = cargaPvDefault(Data, Solares)

    n = size(Data.Red.Branch.T,1);

	Data.Gen.Pv.I = zeros(size(Data.Gen.Tras.I));
	Data.Gen.Pv.pPvg = zeros(size(Data.Gen.Pv.I,1),size(Data.Red.Bus.pCLow,2));
	Data.Gen.Pv.sTop = Data.Gen.Pv.I;
	Data.Gen.Pv.xiTop = Data.Gen.Pv.I;											
	Data.Gen.Pv.pgTop = Data.Gen.Pv.I;
	Data.Gen.Pv.cv = Data.Gen.Pv.I;
	Data.Gen.Pv.cr = Data.Gen.Pv.I;
    
    for i =1:length(Solares)
        [vars] = loadCsvData(Solares(i).fileG,n);

        pPvg = [];
        
        lenVar = length(vars);
        if lenVar == 1
            if strcmp(vars(1).name, 'pPvg') && vars(1).undefBus
                pPvg = vars(1).data;
            end
        end
        if ~isempty(pPvg)
            [sTop_ct, pgTop_ct, cv_ct, cr_ct] = Solares(i).type();
            Data.Gen.Pv.I(Solares(i).nod) = 1;
            Data.Gen.Pv.sTop(Solares(i).nod,:) = sTop_ct;
            Data.Gen.Pv.xiTop(Solares(i).nod,:) = Data.Gen.Pv.sTop(Solares(i).nod,:) .^ 2;											
            Data.Gen.Pv.pgTop(Solares(i).nod,:) = pgTop_ct;
            Data.Gen.Pv.cv(Solares(i).nod,:) = cv_ct;
            Data.Gen.Pv.cr(Solares(i).nod,:) = cr_ct;
            Data.Gen.Pv.pPvg(Solares(i).nod,:) = pPvg;
        end

    end
        

end