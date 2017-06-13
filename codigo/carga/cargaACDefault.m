function [Data] = cargaACDefault(Data, fileTemp, ACs)

    [Data.temp] = cargaTemp(Data, fileTemp);

    n = size(Data.Red.Branch.T,1);
    et = size(Data.Red.Bus.pCLow,2);    

    for i =1:length(ACs)
        [vars] = loadCsvDataSeries(ACs(i).fileT,n);

        tempLow = [];
        tempTop = [];
        tempPref = [];
        beta = [];
        
        lenVar = length(vars);
        if lenVar == 4
            for j = 1:lenVar
                if strcmp(vars(j).name, 'tempLow') && ~vars(j).undefBus
                    tempLow = vars(j).data;
                elseif strcmp(vars(j).name, 'tempTop') && ~vars(j).undefBus
                    tempTop = vars(j).data;
                elseif strcmp(vars(j).name, 'tempPref') && ~vars(j).undefBus
                    tempPref = vars(j).data;
                elseif strcmp(vars(j).name, 'beta') && ~vars(j).undefBus
                    beta = vars(j).data;
                end
            end
        end
        if ~isempty(tempLow) && ~isempty(tempTop) && ~isempty(tempPref) && ~isempty(beta)

            nodsB = max(sign(abs(beta)),[],2);
            nodsTT = max(sign(abs(tempTop)),[],2);
            nodsTL = max(sign(abs(tempLow)),[],2);
            nodsTP = max(sign(abs(tempPref)),[],2);
            
            nods = max(max(max(nodsB,nodsTT),nodsTL),nodsTP);
            
            Data.St.AC.I = nods;
            Data.St.AC.tempLow(:,:) = tempLow;
            Data.St.AC.tempTop(:,:) = tempTop;
            Data.St.AC.tempPref(:,:) = tempPref;
            Data.St.AC.tempIni(:,:) = Data.St.AC.I * ACs(i).tempIni;
            Data.St.AC.epsilon(:,:) = repmat(Data.St.AC.I, [1, et]) * ACs(i).epsilon;
            Data.St.AC.eta(:,:) = repmat(Data.St.AC.I, [1, et]) * ACs(i).eta;
            Data.St.AC.beta(:,:) = beta;
            
            
        end

    end
    
end
