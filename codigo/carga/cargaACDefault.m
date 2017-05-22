function [Data] = cargaACDefault(Data, fileTemp, ACs)

    [Data.temp] = cargaTemp(Data, fileTemp);

    n = size(Data.Red.Branch.T,1);
    et = size(Data.Red.Bus.pCLow,2);    
 
    Data.dt = .25;
    Data.St.AC.I = zeros(size(Data.Red.Branch.T,1),1);
    Data.St.AC.tempLow = zeros(n,et);
    Data.St.AC.tempTop = zeros(n,et);
    Data.St.AC.tempPref = zeros(n,et);
    Data.St.AC.tempIni = Data.St.AC.I;
    Data.St.AC.epsilon = zeros(n,et);
    Data.St.AC.eta = zeros(n,et);
    Data.St.AC.beta = zeros(n,et);


    for i =1:length(ACs)
        [vars] = loadCsvData(ACs(i).fileT,n);

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

            Data.St.AC.I(ACs(1).nod) = 1;
            Data.St.AC.tempLow = tempLow;
            Data.St.AC.tempTop = tempTop;
            Data.St.AC.tempPref = tempPref;
            Data.St.AC.tempIni = Data.St.AC.I * ACs(i).tempIni;
            Data.St.AC.epsilon = repmat(Data.St.AC.I, [1, et]) * ACs(i).epsilon;
            Data.St.AC.eta = repmat(Data.St.AC.I, [1, et]) * ACs(i).eta;
            Data.St.AC.beta = beta;
            
            
        end

    end
    
end
