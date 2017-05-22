function [piPTras, piQmtras, piQMtras] = costosTrasmision(Data,filename)

	n = size(Data.Red.Branch.T,1);
    [vars] = loadCsvData(filename,n);
    lenVar = length(vars);
    piPTras = [];
    piQMtras = [];
    piQmtras = [];
    if lenVar == 3
        for i = 1:lenVar
            if strcmp(vars(i).name, 'piPTras') && vars(i).undefBus
                piPTras = vars(i).data;
            elseif  strcmp(vars(i).name, 'piQMtras') && vars(i).undefBus
                piQMtras = vars(i).data;
            elseif  strcmp(vars(i).name, 'piQmtras') && vars(i).undefBus
                piQmtras = vars(i).data;
            end
        end
        if ~isempty(piPTras) && ~isempty(piQMtras) && ~isempty(piQmtras)
            piPTras = Data.Gen.Tras.I * piPTras;
            piQmtras = Data.Gen.Tras.I * piQmtras;
            piQMtras = Data.Gen.Tras.I * piQMtras;
        end
        
    end
