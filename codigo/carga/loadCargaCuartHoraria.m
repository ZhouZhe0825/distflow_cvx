function [Data] = loadCargaCuartHoraria(Data, filename)

	n = size(Data.Red.Branch.T,1);
    [vars] = loadCsvData(filename,n);
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
        Data.Red.Bus.pCLow = pCLow;
        Data.Red.Bus.qCLow = qCLow;
        ipC = find(sign(sum(abs(Data.Red.Bus.pCLow),2)) == 1);
        iqC = find(sign(sum(abs(Data.Red.Bus.qCLow),2)) == 1);
        Data.Red.Bus.indCons = intersect(ipC, iqC);
    else
        Data.Red.Bus.indCons = [];
    end

% 	n = size(Data.Red.Branch.T,1);
% 	[pCLow, ipC] = loadVarCuartHor(filename, pCSheet, n);
% 	[qCLow, iqC] = loadVarCuartHor(filename, qCSheet, n);
% 
% 	Data.Red.Bus.indCons = intersect(ipC, iqC);
% 
% 	Data.Red.Bus.pCLow = pCLow;
% 	Data.Red.Bus.qCLow = qCLow;
% 
% end
% 
% 
% function [var, indCons] = loadVarCuartHor(filename, sheet, nodos)
% 	[n,~,~] = xlsread(filename, sheet);
% 	var = zeros(nodos,size(n,1)-1);
% 	sirve = ~isnan(n(1,:));
% 	aux = n(:,sirve);
% 	ind = aux(1,:)';
% 	var(ind,:) = aux(2:end,:)';
%     indCons = ind;
% end
% 
