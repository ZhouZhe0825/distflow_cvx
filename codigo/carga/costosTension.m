function [m, delta, cdv] = costosTension(Data,filename)

	n = size(Data.Red.Branch.T,1);
    [vars] = loadCsvData(filename,n);
    lenVar = length(vars);
    m = [];
    delta = [];
    cdv = [];
    if lenVar == 3
        for i = 1:lenVar
            if strcmp(vars(i).name, 'm') && ~vars(i).undefBus
                m = vars(i).data;
            elseif strcmp(vars(i).name, 'delta') && ~vars(i).undefBus
                delta = vars(i).data;
            elseif strcmp(vars(i).name, 'cdv') && ~vars(i).undefBus
                cdv = vars(i).data;
            end
        end
    end
