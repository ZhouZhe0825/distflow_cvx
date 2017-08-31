function [Data] = cargaTrasDefault(Data, Tras)

    n = size(Data.Red.Branch.T,1);

    for i =1:length(Tras)
        [vars] = loadCsvDataSeries(Tras(i).fileG,n);

        pgLow = [];
        pgTop = [];
        qgLow = [];
        qgTop = [];
        
        lenVar = length(vars);
        if lenVar == 4
            for j = 1:lenVar
                if strcmp(vars(j).name, 'pgLow') && vars(j).undefBus
                    pgLow = vars(j).data;
                elseif strcmp(vars(j).name, 'pgTop') && vars(j).undefBus
                    pgTop = vars(j).data;
                elseif strcmp(vars(j).name, 'qgLow') && vars(j).undefBus
                    qgLow = vars(j).data;
                elseif strcmp(vars(j).name, 'qgTop') && vars(j).undefBus
                    qgTop = vars(j).data;
                end
            end
        end
        if ~isempty(pgLow) && ~isempty(qgLow) && ~isempty(pgTop) && ~isempty(qgTop)
            Data.Red.Bus.uLow(Tras(i).nod,:) = Tras(i).uLow;
            Data.Red.Bus.uTop(Tras(i).nod,:) = Tras(i).uTop;
            Data.Gen.Tras.pgLow(Tras(i).nod,:) = pgLow;
            Data.Gen.Tras.qgLow(Tras(i).nod,:) = qgLow;
            Data.Gen.Tras.pgTop(Tras(i).nod,:) = pgTop;
            Data.Gen.Tras.qgTop(Tras(i).nod,:) = qgTop;
            Data.Gen.Tras.I(Tras(i).nod,:) = 1;
        end

    end
        

end