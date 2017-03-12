function [Data] = loadCargaCuartHorariaApp(filename, Data, App, pCSheet, qCSheet, minpC, minqC)

	n = size(Data.Red.Branch.T,1);
    
    [aux, ~] = loadVarCuartHor(filename, [pCSheet, '_1'], n, minpC);

    Data.Util.pzCnPref = zeros(size(aux,1), size(aux,2), length(App));
    Data.Util.pzCnLow = Data.Util.pzCnPref;
    Data.Util.pzCnTop = Data.Util.pzCnPref;
    Data.Red.Bus.indCons = [];
    
    pCLow = Data.Util.pzCnPref;
    qCLow = pCLow;
    
    for a = 1:length(App)
        app = App(a);
        
        pCSheetApp = [pCSheet, '_', num2str(a)];
        qCSheetApp = [qCSheet, '_', num2str(a)];
        
        [pCLowApp, ipCApp] = loadVarCuartHor(filename, pCSheetApp, n, minpC);
        [qCLowApp, iqCApp] = loadVarCuartHor(filename, qCSheetApp, n, minqC);
        
        Data.Util.pzCnPref(:,:,a) = pCLowApp;
        Data.Util.pzCnLow(:,:,a) = pCLowApp * app.nMultipLow;
        Data.Util.pzCnTop(:,:,a) = pCLowApp * app.nMultipTop;
        
        pCLow = pCLow + pCLowApp;
        qCLow = qCLow + qCLowApp;
    	Data.Red.Bus.indCons = intersect(intersect(ipCApp, iqCApp), Data.Red.Bus.indCons);
    end
    
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

