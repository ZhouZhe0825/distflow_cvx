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
        [vars] = loadCsvDataSeries(Solares(i).fileG,n);

        pPvg = [];
        
        lenVar = length(vars);
        if lenVar == 1
            if strcmp(vars(1).name, 'pPvg') && vars(1).undefBus
                pPvg = vars(1).data;
            end
        end
        if ~isempty(pPvg)
            fields = {'sTop_ct', 'pgTop_ct', 'cv_ct', 'cr_ct'};
            S = loadCsvDef(Solares(i).type(), fields);
            Data.Gen.Pv.I(Solares(i).nod) = 1;
            Data.Gen.Pv.sTop(Solares(i).nod,:) = S.sTop_ct;
            Data.Gen.Pv.xiTop(Solares(i).nod,:) = Data.Gen.Pv.sTop(Solares(i).nod,:) .^ 2;											
            Data.Gen.Pv.pgTop(Solares(i).nod,:) = S.pgTop_ct;
            Data.Gen.Pv.cv(Solares(i).nod,:) = S.cv_ct;
            Data.Gen.Pv.cr(Solares(i).nod,:) = S.cr_ct;
            Data.Gen.Pv.pPvg(Solares(i).nod,:) = pPvg;

            
        end

    end
        

end