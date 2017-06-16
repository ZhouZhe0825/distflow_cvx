function [Header] = createHeader(Var, Data, Config)

	horas = loadHoras();

    [rowTap, colTap, ~] = find(triu(Data.Red.Branch.Itap) == 1);
    cantTaps = length(rowTap);

    indCaps = find(Data.Red.Bus.Icap == 1);
    cantCaps = length(indCaps);

    [rowTReg, colTReg, ~] = find(triu(Data.Red.Branch.Itreg) == 1);
    cantReg = length(rowTReg);

	nodCh = find(Data.ClNI.I == 1);
    cantCarg = length(nodCh);
	indHeadEt = (2:1+Config.Etapas);

	n = size(Data.Red.Branch.T,1);
	nodos = (1:n)';
    
	Header.Main = cell(1+cantTaps+cantCaps+cantCarg,Config.Etapas+1);

	Header.Main(1,indHeadEt) = horas(Config.iniEtapa:Config.iniEtapa+Config.Etapas-1);

    if cantTaps > 0
        for i = 1:cantTaps
            Header.Main{1+i,1} = ['Ntr_' num2str(rowTap(i)) '-' num2str(colTap(i))];
            Header.Main(1+i,indHeadEt) = num2cell(squeeze(round(Var.Red.Branch.Ntr(rowTap(i),colTap(i),:)))');
        end
    end

    if cantReg > 0
        for i = 1:cantReg
            Header.Main{1+cantTaps+i,1} = ['Reg_' num2str(rowTReg(i)) '-' num2str(colTReg(i))];
            Header.Main(1+cantTaps+i,indHeadEt) = num2cell(squeeze(Var.Red.Branch.Rtr(rowTReg(i),colTReg(i),:))');
        end
    end

    if cantCaps > 0
        for i = 1:cantCaps
            Header.Main{1+cantTaps+cantReg+i,1} = ['Ncp_n_' num2str(indCaps(i))];
        end
        Header.Main((2+cantReg+cantTaps:1+cantTaps+cantReg+cantCaps),indHeadEt) = num2cell(squeeze(round(Var.Red.Bus.Ncp(indCaps,:,:)))');
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


