function [c] = costosGBasico(Data,GenBas)

	n = size(Data.Red.Branch.T,1);
    et = size(Data.Red.Bus.pCLow,2);
    c = zeros(n,et);
    for i = 1:length(GenBas)
        [vars] = loadCsvDataSeries(GenBas(i).fileC,n);
        lenVar = length(vars);
        c_aux = [];
        if lenVar == 3
            for j = 1:lenVar
                if strcmp(vars(j).name, 'c') && vars(j).undefBus
                    c_aux = vars(j).data;
                end
            end
        end
        if ~isempty(c_aux)
            c(GenBas(i).nod,:) = c_aux;
        end
        
    end
