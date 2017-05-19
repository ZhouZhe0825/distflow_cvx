function [Data] = cargaUtilDefault(Data, tgPhi, betaE, betaT, Cargas, App)

Data.Util.tgPhi = tgPhi;

Data.Util.betaT = zeros(size(Data.Red.Bus.pCLow,1),size(betaT,1));
Data.Util.betaT(Data.Red.Bus.indCons,:) = ones(length(Data.Red.Bus.indCons),1)*betaT';

Data.Util.betaE = zeros(size(Data.Red.Bus.pCLow));
for estCg = 1:size(Cargas,1)
	Data.Util.betaE(Cargas(estCg).nod,:) = betaE;
end

Data.Util.pzCnPref = repmat(full(Data.Red.Bus.pCLow), [1,1,length(App)]);
Data.Util.pzCnLow = Data.Util.pzCnPref;
Data.Util.pzCnTop = Data.Util.pzCnPref;

for a = 1:length(App)
    app = App(a);
    Data.Util.pzCnPref(:,:,app.I) = Data.Util.pzCnPref(:,:,app.I) * app.Pref;
    Data.Util.pzCnLow(:,:,app.I) = Data.Util.pzCnLow(:,:,app.I) * app.Low * app.nMultipLow;
    Data.Util.pzCnTop(:,:,app.I) = Data.Util.pzCnTop(:,:,app.I) * app.Top * app.nMultipTop;
end

Data.Util.pzCnPrefE = Data.ClNI.pC;
Data.Util.pzCnLowE = Data.Util.pzCnPrefE .* Data.ClNI.nMultipLow;
Data.Util.pzCnTopE = Data.Util.pzCnPrefE .* Data.ClNI.nMultipTop;

Data.Util.qzCnPrefE = Data.ClNI.qC;
Data.Util.qzCnLowE = Data.Util.qzCnPrefE .* Data.ClNI.nMultipLow;
Data.Util.qzCnTopE = Data.Util.qzCnPrefE .* Data.ClNI.nMultipTop;
