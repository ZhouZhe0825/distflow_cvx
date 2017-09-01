function [Data] = cargaGBasicoDefault(Data, GenBas)

    n = size(Data.Red.Branch.T,1);

    for i =1:length(GenBas)
        [vars] = loadCsvDataSeries(GenBas(i).fileG,n);

        pgLow = [];
        qgLow = [];
        pgTop = [];
        qgTop = [];
        
        lenVar = length(vars);
        if lenVar == 4
            for j = 1:lenVar
                if strcmp(vars(j).name, 'pgLow') && vars(j).undefBus
                    pgLow = vars(j).data;
                elseif strcmp(vars(j).name, 'qgLow') && vars(j).undefBus
                    qgLow = vars(j).data;
                elseif strcmp(vars(j).name, 'pgTop') && vars(j).undefBus
                    pgTop = vars(j).data;
                elseif strcmp(vars(j).name, 'qgTop') && vars(j).undefBus
                    qgTop = vars(j).data;
                end
            end
        end
        if ~isempty(pgLow) && ~isempty(qgLow) && ~isempty(pgTop) && ~isempty(qgTop)
            Data.Gen.Basic.pgLow(GenBas(i).nod,:) = pgLow;
            Data.Gen.Basic.qgLow(GenBas(i).nod,:) = qgLow;
            Data.Gen.Basic.pgTop(GenBas(i).nod,:) = pgTop;
            Data.Gen.Basic.qgTop(GenBas(i).nod,:) = qgTop;
            Data.Gen.Basic.pgIni(GenBas(i).nod,:) = GenBas(i).pgIni;
            Data.Gen.Basic.qgIni(GenBas(i).nod,:) = GenBas(i).qgIni;
            Data.Gen.Basic.I(GenBas(i).nod,:) = 1;
        end

    end
        

end