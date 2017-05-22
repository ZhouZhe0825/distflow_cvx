function [Data] = cargaBatDefault(Data, Baterias)

    n = size(Data.Red.Branch.T,1);

    
	Data.St.Bat.I = zeros(size(Data.Red.Branch.T,1),1);
	Data.St.Bat.cv = Data.St.Bat.I;
	Data.St.Bat.cr = Data.St.Bat.I;
	Data.St.Bat.epsilon = Data.St.Bat.I;
	Data.St.Bat.eta = Data.St.Bat.I;
	Data.St.Bat.pgTop = Data.St.Bat.I;
	Data.St.Bat.pgLow = Data.St.Bat.I;
	Data.St.Bat.sTop = Data.St.Bat.I;
	Data.St.Bat.ETop = Data.St.Bat.I;
	Data.St.Bat.ELow = Data.St.Bat.I;
	Data.St.Bat.beta = Data.St.Bat.I;
	Data.St.Bat.wU = Data.St.Bat.I;
	Data.St.Bat.wOm = Data.St.Bat.I;
	Data.St.Bat.m1 = Data.St.Bat.I;
	Data.St.Bat.m2 = Data.St.Bat.I;
	Data.St.Bat.m3 = Data.St.Bat.I;
    Data.St.Bat.EIni = Data.St.Bat.I;
    Data.St.Bat.xiTop = Data.St.Bat.I;
    Data.St.Bat.kapa = Data.St.Bat.I;
    Data.St.Bat.gama = Data.St.Bat.I;

    for i =1:length(Baterias)
            [cv,cr,epsilon,eta,pgTop,pgLow,sTop,ETop,ELow,beta,wU,wOm,m1,m2,m3,kapa,gama] = Baterias(i).type(); 
        
            Data.St.Bat.I(Baterias(i).nod) = 1;
            Data.St.Bat.cv(Baterias(i).nod) = cv;
            Data.St.Bat.cr(Baterias(i).nod) = cr;
            Data.St.Bat.epsilon(Baterias(i).nod) = epsilon;
            Data.St.Bat.eta(Baterias(i).nod) = eta;
            Data.St.Bat.pgTop(Baterias(i).nod) = pgTop;
            Data.St.Bat.pgLow(Baterias(i).nod) = pgLow;
            Data.St.Bat.sTop(Baterias(i).nod) = sTop;
            Data.St.Bat.ETop(Baterias(i).nod) = ETop;
            Data.St.Bat.ELow(Baterias(i).nod) = ELow;
            Data.St.Bat.beta(Baterias(i).nod) = beta;
            Data.St.Bat.wU(Baterias(i).nod) = wU;
            Data.St.Bat.wOm(Baterias(i).nod) = wOm;
            Data.St.Bat.m1(Baterias(i).nod) = m1;
            Data.St.Bat.m2(Baterias(i).nod) = m2;
            Data.St.Bat.m3(Baterias(i).nod) = m3;
            Data.St.Bat.EIni(Baterias(i).nod) = Baterias(i).EIni;
            Data.St.Bat.xiTop(Baterias(i).nod) = Data.St.Bat.sTop(Baterias(i).nod).^2;
            Data.St.Bat.kapa(Baterias(i).nod) = kapa;
            Data.St.Bat.gama(Baterias(i).nod) = gama;
    end
        

end