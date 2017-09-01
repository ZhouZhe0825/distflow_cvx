function [Data] = cargaPvDefault(Data, Solares)

    n = size(Data.Red.Branch.T,1);

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
            Data.Gen.Pv.pgIni(Solares(i).nod,:) = Solares(i).pgIni;
            Data.Gen.Pv.qgIni(Solares(i).nod,:) = Solares(i).qgIni;

            
        end

    end
        

end