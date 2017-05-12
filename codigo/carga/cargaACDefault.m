function [Data] = cargaACDefault(Data)

    Data.dt = .25;
    Data.St.AC.I = zeros(size(Data.Red.Branch.T,1),1);
    Data.St.AC.tempLow = repmat(Data.St.AC.I, [1, size(Data.temp,2)]);
    Data.St.AC.tempTop = repmat(Data.St.AC.I, [1, size(Data.temp,2)]);
    Data.St.AC.tempPref = repmat(Data.St.AC.I, [1, size(Data.temp,2)]);
    Data.St.AC.tempIni = Data.St.AC.I;
    Data.St.AC.epsilon = repmat(Data.St.AC.I, [1, size(Data.temp,2)]);
    Data.St.AC.eta = repmat(Data.St.AC.I, [1, size(Data.temp,2)]);
    Data.St.AC.beta = repmat(Data.St.AC.I, [1, size(Data.temp,2)]);
end
