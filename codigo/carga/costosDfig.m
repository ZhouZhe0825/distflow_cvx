function [rhopWi,rhomqWi,rhoMqWi] = costosDfig(Data,Eolicos)

	n = size(Data.Red.Branch.T,1);
    et = size(Data.Red.Bus.pCLow,2);
    rhopWi = zeros(n,et);
    rhomqWi = zeros(n,et);
    rhoMqWi = zeros(n,et);
    for i = 1:length(Eolicos)
        [vars] = loadCsvDataSeries(Eolicos(i).fileC,n);
        lenVar = length(vars);
        rhopWi_aux = [];
        rhomqWi_aux = [];
        rhoMqWi_aux = [];
        if lenVar == 3
            for j = 1:lenVar
                if strcmp(vars(j).name, 'rhopWi') && vars(j).undefBus
                    rhopWi_aux = vars(j).data;
                elseif strcmp(vars(j).name, 'rhomqWi') && vars(j).undefBus
                    rhomqWi_aux = vars(j).data;
                elseif strcmp(vars(j).name, 'rhoMqWi') && vars(j).undefBus
                    rhoMqWi_aux = vars(j).data;
                end
            end
        end
        if ~isempty(rhopWi_aux) && ~isempty(rhomqWi_aux) && ~isempty(rhoMqWi_aux)
            rhopWi(Eolicos(i).nod,:) = rhopWi_aux;
            rhomqWi(Eolicos(i).nod,:) = rhomqWi_aux;											
            rhoMqWi(Eolicos(i).nod,:) = rhoMqWi_aux;
        end
        
    end
