function [betaE] = utilBetasE(Data,Cargas)

	n = size(Data.Red.Branch.T,1);
    et = size(Data.Red.Bus.pCLow,2);
    betaE = zeros(n,et);
    for i = 1:length(Cargas)
        [vars] = loadCsvDataSeries(Cargas(i).fileU,n);
        lenVar = length(vars);
        betaE_aux = [];
        if lenVar == 1
            if strcmp(vars(1).name, 'betaE') && vars(1).undefBus
                betaE_aux = vars(1).data;
            end
        end
        if ~isempty(betaE_aux)
            betaE(Cargas(i).nod,:) = betaE_aux;
        end
    end
