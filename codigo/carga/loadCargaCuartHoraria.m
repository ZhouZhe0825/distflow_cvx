function [Data] = loadCargaCuartHoraria(filename, Data, pCSheet, qCSheet, minpC, minqC)

	n = size(Data.Red.Branch.T,1);
	[pCLow, ipC] = loadVarCuartHor(filename, pCSheet, n, minpC);
	[qCLow, iqC] = loadVarCuartHor(filename, qCSheet, n, minqC);

	Data.Red.Bus.indCons = intersect(ipC, iqC);

	Data.Red.Bus.pCLow = pCLow;
	Data.Red.Bus.qCLow = qCLow;

end


function [var, indCons] = loadVarCuartHor(filename, sheet, nodos, minVal)
	[n,t,r] = xlsread(filename, sheet);
	var = ones(nodos,size(n,1)-1)*minVal;
	sirve = ~isnan(n(1,:));
	aux = n(:,sirve);
	ind = aux(1,:)';
	var(ind,:) = aux(2:end,:)';
    indCons = ind;
end

