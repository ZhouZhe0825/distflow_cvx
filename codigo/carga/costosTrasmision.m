function [piPTras, piQmtras, piQMtras] = costosTrasmision(Data,Tras)

	n = size(Data.Red.Branch.T,1);
    et = size(Data.Red.Bus.pCLow,2);
    piPTras = zeros(n,et);
    piQMtras = zeros(n,et);
    piQmtras = zeros(n,et);
    for i = 1:length(Tras)
        [vars] = loadCsvDataSeries(Tras(i).fileC,n);
        lenVar = length(vars);
        piPTras_aux = [];
        piQMtras_aux = [];
        piQmtras_aux = [];
        if lenVar == 3
            for j = 1:lenVar
                if strcmp(vars(j).name, 'piPTras') && vars(j).undefBus
                    piPTras_aux = vars(j).data;
                elseif strcmp(vars(j).name, 'piQMtras') && vars(j).undefBus
                    piQMtras_aux = vars(j).data;
                elseif strcmp(vars(j).name, 'piQmtras') && vars(j).undefBus
                    piQmtras_aux = vars(j).data;
                end
            end
        end
        if ~isempty(piPTras_aux) && ~isempty(piQMtras_aux) && ~isempty(piQmtras_aux)
            piPTras(Tras(i).nod,:) = piPTras_aux;
            piQMtras(Tras(i).nod,:) = piQMtras_aux;											
            piQmtras(Tras(i).nod,:) = piQmtras_aux;
        end
        
    end
    