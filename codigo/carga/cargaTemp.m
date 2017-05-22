function [temp] = cargaTemp(Data, filename)

    n = size(Data.Red.Branch.T,1);
    temp = [];
    [vars] = loadCsvData(filename,n);
    lenVar = length(vars);
    if lenVar == 1
        if strcmp(vars(1).name, 'temp') && vars(1).undefBus
            temp = vars(1).data;
        end
    end
