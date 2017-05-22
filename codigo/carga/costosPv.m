function [rhopPv,rhomqPv,rhoMqPv] = costosPv(Data,Solares)

	n = size(Data.Red.Branch.T,1);
    et = size(Data.Red.Bus.pCLow,2);
    rhopPv = zeros(n,et);
    rhomqPv = zeros(n,et);
    rhoMqPv = zeros(n,et);
    for i = 1:length(Solares)
        [vars] = loadCsvData(Solares(i).fileC,n);
        lenVar = length(vars);
        rhopPv_aux = [];
        rhomqPv_aux = [];
        rhoMqPv_aux = [];
        if lenVar == 3
            for j = 1:lenVar
                if strcmp(vars(j).name, 'rhopPv') && vars(j).undefBus
                    rhopPv_aux = vars(j).data;
                elseif strcmp(vars(j).name, 'rhomqPv') && vars(j).undefBus
                    rhomqPv_aux = vars(j).data;
                elseif strcmp(vars(j).name, 'rhoMqPv') && vars(j).undefBus
                    rhoMqPv_aux = vars(j).data;
                end
            end
        end
        if ~isempty(rhopPv_aux) && ~isempty(rhomqPv_aux) && ~isempty(rhoMqPv_aux)
            rhopPv(Solares(i).nod,:) = rhopPv_aux;
            rhomqPv(Solares(i).nod,:) = rhomqPv_aux;											
            rhoMqPv(Solares(i).nod,:) = rhoMqPv_aux;
        end
        
    end
