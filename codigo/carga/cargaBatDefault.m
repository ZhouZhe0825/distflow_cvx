function [Data] = cargaBatDefault(Data, Baterias)

    for i =1:length(Baterias)
            fields = {...
                'cv','cr','epsilon','eta','pgTop','pgLow','sTop','ETop',...
                'ELow','beta','wU','wOm','m1','m2','m3','kapa','gama'...
                };
            S = loadCsvDef(Baterias(i).type(), fields);
        
            Data.St.Bat.I(Baterias(i).nod) = 1;
            Data.St.Bat.cv(Baterias(i).nod) = S.cv;
            Data.St.Bat.cr(Baterias(i).nod) = S.cr;
            Data.St.Bat.epsilon(Baterias(i).nod) = S.epsilon;
            Data.St.Bat.eta(Baterias(i).nod) = S.eta;
            Data.St.Bat.pgTop(Baterias(i).nod) = S.pgTop;
            Data.St.Bat.pgLow(Baterias(i).nod) = S.pgLow;
            Data.St.Bat.sTop(Baterias(i).nod) = S.sTop;
            Data.St.Bat.ETop(Baterias(i).nod) = S.ETop;
            Data.St.Bat.ELow(Baterias(i).nod) = S.ELow;
            Data.St.Bat.beta(Baterias(i).nod) = S.beta;
            Data.St.Bat.wU(Baterias(i).nod) = S.wU;
            Data.St.Bat.wOm(Baterias(i).nod) = S.wOm;
            Data.St.Bat.m1(Baterias(i).nod) = S.m1;
            Data.St.Bat.m2(Baterias(i).nod) = S.m2;
            Data.St.Bat.m3(Baterias(i).nod) = S.m3;
            Data.St.Bat.EIni(Baterias(i).nod) = Baterias(i).EIni;
            Data.St.Bat.xiTop(Baterias(i).nod) = Data.St.Bat.sTop(Baterias(i).nod).^2;
            Data.St.Bat.kapa(Baterias(i).nod) = S.kapa;
            Data.St.Bat.gama(Baterias(i).nod) = S.gama;
    end
        

end