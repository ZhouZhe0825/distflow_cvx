function [Data] = loadCargaCuartHoraria(filename, Data, pCSheet, qCSheet)

	n = size(Data.Red.Branch.T,1);
	[pCLow, ipC] = loadVarCuartHor(filename, pCSheet, n);
	[qCLow, iqC] = loadVarCuartHor(filename, qCSheet, n);

	Data.Red.Bus.indCons = intersect(ipC, iqC);

	Data.Red.Bus.pCLow = pCLow;
	Data.Red.Bus.qCLow = qCLow;

end


function [var, indCons] = loadVarCuartHor(filename, sheet, nodos)
	[n,~,~] = xlsread(filename, sheet);
	var = zeros(nodos,size(n,1)-1);
	sirve = ~isnan(n(1,:));
	aux = n(:,sirve);
	ind = aux(1,:)';
	var(ind,:) = aux(2:end,:)';
    indCons = ind;
end

