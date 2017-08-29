function [Data] = cargaUtilDefault(Data, fileUtilBetaT, utilOptFuncCuad, Cargas, App)

[Data.Util.betaE] = utilBetasE(Data,Cargas);

[Data.Util.betaT] = utilBetasT(Data,fileUtilBetaT);

Data.Util.betaTcuad = utilOptFuncCuad;

Data.Util.pzCnPref = repmat(full(Data.Red.Bus.pCLow), [1,1,length(App)]);
 
Data.Util.pzCnLow = Data.Util.pzCnPref;
 
Data.Util.pzCnTop = Data.Util.pzCnPref;
 
for a = 1:length(App)
    app = App(a);
    Data.Util.pzCnPref(:,:,app.I) = Data.Red.Bus.pCLow(:,:) * app.Pref;
    Data.Util.pzCnLow(:,:,app.I) = Data.Red.Bus.pCLow(:,:) * app.Low * app.nMultipLow;
    Data.Util.pzCnTop(:,:,app.I) = Data.Red.Bus.pCLow(:,:) * app.Top * app.nMultipTop;
    Data.Util.tgPhi(a) = app.tgPhi;
end

Data.Util.pzCnPrefE(:,:) = Data.ClNI.pC;
Data.Util.pzCnLowE(:,:) = Data.Util.pzCnPrefE .* Data.ClNI.nMultipLow;
Data.Util.pzCnTopE(:,:) = Data.Util.pzCnPrefE .* Data.ClNI.nMultipTop;

Data.Util.qzCnPrefE(:,:) = Data.ClNI.qC;
Data.Util.qzCnLowE(:,:) = Data.Util.qzCnPrefE .* Data.ClNI.nMultipLow;
Data.Util.qzCnTopE(:,:) = Data.Util.qzCnPrefE .* Data.ClNI.nMultipTop;
