function [vars] = loadCsvDataSeries(filename,n)

aux = importdata(filename);
csv.nums = aux.data;
csv.header = aux.textdata(1,2:end);
csv.hora = aux.textdata(2:end,1);
splheader = cell(length(csv.header),2);
for i = 1:length(csv.header)
	aux = strsplit(csv.header{i},'@');
	splheader{i,1} = aux{1};
	splheader{i,2} = aux{2};	
end
uniqueVars = unique(splheader(:,1));
vars = [];
for i = 1:length(uniqueVars)
	var.name = uniqueVars{i};
	var.data = zeros(n,length(csv.hora));
    var.undefBus = true;
	vars = [vars;var];
end
for i = 1:size(splheader,1)
	j = 1;
	while j <= length(vars) && strcmp(splheader{i,1},vars(j).name) == 0
		j = j+1;
	end
	if strcmp(splheader{i,1},vars(j).name)
		textInd = splheader{i,2};
		if isempty(textInd) % para nodo a definir
            vars(j).data = csv.nums(:,i)';
		elseif strcmp(textInd,'all') % para todos los nodos
            vars(j).data = ones(n,1)*csv.nums(:,i)';
            vars(j).undefBus = false;
		else
			[ind, status] = str2num(textInd);
			if status && ind <= n % para indice
				vars(j).data(ind,:) = csv.nums(:,i)';
                vars(j).undefBus = false;
			end
		end
	end
end