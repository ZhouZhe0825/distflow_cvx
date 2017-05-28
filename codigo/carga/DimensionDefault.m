function [n,et,app] = DimensionDefault(inFilename,bus_sheet,fileCurvaCarga,App)

    [n_bu,~,~] = xlsread(inFilename, bus_sheet);
    
    aux = importdata(fileCurvaCarga);


    n = length(n_bu(:,1));
    et = size(aux.data,1);
    app = size(App,1);

end
	