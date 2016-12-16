function [pCLowInv, pCLowVer, qCLowInv, qCLowVer, iCi, iCv] = loadCargaCuartHoraria(filename, nodos, minpC, minqC)
% filename = 'carga_subredLP.xlsx';
pCLowShInv = 'Ppu_inv';
pCLowShVer = 'Ppu_ver';
qCLowShInv = 'Qpu_inv';
qCLowShVer = 'Qpu_ver';


[pCLowInv, ipCi] = loadVarCuartHor(filename, pCLowShInv, nodos, minpC);
[pCLowVer, iqCi] = loadVarCuartHor(filename, pCLowShVer, nodos, minpC);
[qCLowInv, ipCv] = loadVarCuartHor(filename, qCLowShInv, nodos, minqC);
[qCLowVer, iqCv] = loadVarCuartHor(filename, qCLowShVer, nodos, minqC);

iCi = intersect(ipCi, iqCi);
iCv = intersect(ipCv, iqCv);

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

