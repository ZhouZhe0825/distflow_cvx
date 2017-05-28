function [Data] = loadCargaCuartHoraria(Data, filename)

	n = size(Data.Red.Branch.T,1);
    [vars] = loadCsvDataSeries(filename,n);
    lenVar = length(vars);
    pCLow = [];
    qCLow = [];
    if length(vars) == 2
        for i = 1:lenVar
            if strcmp(vars(i).name, 'pCLow') && ~vars(i).undefBus
                pCLow = vars(i).data;
            elseif strcmp(vars(i).name, 'qCLow') && ~vars(i).undefBus
                qCLow = vars(i).data;
            end
        end
    end
    if ~isempty(pCLow) && ~isempty(qCLow)
        Data.Red.Bus.pCLow(:,:) = pCLow;
        Data.Red.Bus.qCLow(:,:) = qCLow;
        ipC = find(sign(sum(abs(Data.Red.Bus.pCLow),2)) == 1);
        iqC = find(sign(sum(abs(Data.Red.Bus.qCLow),2)) == 1);
        Data.Red.Bus.Icons(intersect(ipC, iqC)) = 1;
    end
end