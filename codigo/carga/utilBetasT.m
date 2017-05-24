function [betaT] = utilBetasT(Data,filename)

	n = size(Data.Red.Branch.T,1);
    et = size(Data.Red.Bus.pCLow,2);
    betaT = zeros(n,et);
    [vars] = loadCsvDataSeries(filename,n);
    lenVar = length(vars);
    betaT_aux = [];
    if lenVar == 1
        if strcmp(vars(1).name, 'betaT') && vars(1).undefBus
            betaT_aux = vars(1).data;
        end
    end
    if ~isempty(betaT_aux)
        betaT(Data.Red.Bus.indCons,:) = ones(size(Data.Red.Bus.indCons))*betaT_aux;
    end
