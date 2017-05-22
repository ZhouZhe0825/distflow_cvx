function [Header] = createHeader(Var, Data, Config, cantTaps, cantCaps, cantCarg)

	horas = loadHoras();

	Taps = matOverTime(abs(Data.Red.Bus.Tap + Data.Red.Bus.NtrLow + Data.Red.Bus.NtrTop));
	indTaps = find(Taps == 1);
	Caps = matOverTime(abs(Data.Red.Bus.Ncp + Data.Red.Bus.CapLow + Data.Red.Bus.CapTop));
	indCaps = find(Caps == 1);

	nodCh = find(Data.ClNI.I == 1);
	indHeadEt = (2:1+Config.Etapas);

	n = size(Data.Red.Branch.T,1);
	nodos = (1:n)';
    
	Header.Main = cell(1+cantTaps+cantCaps+cantCarg,Config.Etapas+1);

	Header.Main(1,indHeadEt) = horas(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1);

    if cantTaps > 0
        for i = 1:cantTaps
            Header.Main{1+i,1} = ['Ntr' num2str(i) '_n_' num2str(indTaps(i))];
        end
        Header.Main((2:1+cantTaps),indHeadEt) = num2cell(squeeze(round(Var.Red.Bus.Ntr(indTaps,:,:)))');
    end

    if cantCaps > 0
        for i = 1:cantCaps
            Header.Main{1+cantTaps+i,1} = ['Cap' num2str(i) '_n_' num2str(indCaps(i))];
        end
        Header.Main((2+cantTaps:1+cantTaps+cantCaps),indHeadEt) = num2cell(squeeze(round(Var.Red.Bus.Cap(indCaps,:,:)))');
    end

    if cantCarg > 0
        for i = 1:cantCarg
            Header.Main{1+cantTaps+cantCaps+i,1} = ['Chg_' num2str(Data.ClNI.pC(nodCh(i))) '_d_' num2str(Data.ClNI.d(nodCh(i))) '_n_' num2str(nodCh(i))];
        end
        Header.Main((2+cantTaps+cantCaps:1+cantTaps+cantCaps+cantCarg),indHeadEt) = num2cell(squeeze(round(Var.ClNI.on(nodCh,:))));
    end

	TotalT = matOverTime(Data.Red.Branch.T);
    [rowT, colT, ~] = find(TotalT == 1);
	
	Header.Bus = strrep(mat2cell([repmat(['bus_'],length(nodos),1), num2str(nodos)], ones(1,length(nodos)), 4+(ceil(log(n)/log(10)))), ' ', '');
    Header.Branch = strrep(mat2cell([repmat(['br_'],length(rowT),1), num2str(rowT), repmat(['-'],length(rowT),1), num2str(colT)], ones(1,length(rowT)), 4+2*(ceil(log(n)/log(10)))), ' ', '');
    
end


