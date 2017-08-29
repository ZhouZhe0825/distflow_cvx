function printCost_(Header, Var, Data, outFilename)

    et = size(Var.Red.Bus.v,2);

    cV       = zeros(1,et);
    cPtr     = zeros(1,et);
    cQtr     = zeros(1,et);
    cPgen    = zeros(1,et);
    cQgen    = zeros(1,et);
    cF       = zeros(1,et);
    cTr      = zeros(1,et);
    cSw      = zeros(1,et);
    cCp      = zeros(1,et);
    cBatStb  = zeros(1,et);
    cBatBeta = zeros(1,et);
	cClNI    = zeros(1,et);
	cPpv     = zeros(1,et);
	cQpv     = zeros(1,et);
	cPwi     = zeros(1,et);
	cQwi     = zeros(1,et);
	cAC      = zeros(1,et);

    cV =       sum(Data.Cost.cdv.*Var.Red.Bus.cDv,1);
    cPtr =     sum(Data.Cost.piPTras.*Var.Red.Bus.PTras,1);
    cQtr =     sum(Var.Red.Bus.cQTras,1);
    
    if Data.Util.betaTcuad
        cF = sum(Data.Util.betaT(:,:,1).*((Var.Red.Bus.pCn(:,:,1) - Data.Util.pzCnPref(:,:,1)).^2),1);
    else
        cF = sum(Data.Util.betaT(:,:,1).*(Data.Util.pzCnPref(:,:,1)-(Var.Red.Bus.pCn(:,:,1))),1);
    end

    cTr =      sum(Data.Cost.cTap.*(Var.Red.Branch.NtrDif.^2),1);
    cSw =      sum(Data.Cost.cY.*(Var.Red.Branch.yDif.^2),1);
    cCp =      sum(Data.Cost.cCap.*(Var.Red.Bus.NcpDif.^2),1);

	if isfield(Var, 'ClNI')
		cClNI =    sum(Data.Util.betaE.*((Var.ClNI.pC - Data.Util.pzCnPrefE).^2),1);
    end
    

	if isfield(Var, 'ClRes')
        if isfield(Var.ClRes, 'Tvar')
    		cAC =    sum(Data.St.AC.beta.*(Var.ClRes.Tvar - Data.St.AC.tempPref).^2,1);
        end
	end

	if isfield(Var, 'Gen')
		if isfield(Var.Gen, 'Basic')
            cPgen =    sum(Data.Cost.cBas .* Var.Gen.Basic.pGBas,1);
            cQgen =    sum(Var.Gen.Basic.cqGBas);
        end
		if isfield(Var.Gen, 'Pv')
			cPpv =     sum(Data.Cost.rhopPv.*Var.Gen.Pv.pPv,1);
			cQpv =     sum(Var.Gen.Pv.cqPv,1);
		end
		if isfield(Var.Gen, 'Dfig')
			cPwi =    sum(Data.Cost.rhopWi .* Var.Gen.Dfig.pWi,1);
			cQwi =    sum(Var.Gen.Dfig.cqWi,1);
		end
    end
	if isfield(Var, 'St')
		if isfield(Var.St, 'Bat')
            cBatStb =  sum(Var.St.Bat.cStb,1);
            cBatBeta = sum(Data.St.Bat.beta.* ((Data.St.Bat.EPref - Var.St.Bat.EStb).^2) + Data.St.Bat.wU,1);
        end
    end

    rowHeader = cell(17,1);
    rowHeader{1} =  'cV';
    rowHeader{2} =  'cPtr';
    rowHeader{3} =  'cQtr';
    rowHeader{4} =  'cPgen';
    rowHeader{5} =  'cQgen';
    rowHeader{6} =  'cF';
    rowHeader{7} =  'cTr';
    rowHeader{8} =  'cSw';
    rowHeader{9} =  'cCp';
    rowHeader{10} = 'cBatStb';
    rowHeader{11} = 'cBatBeta';
    rowHeader{12} = 'cPpv';
    rowHeader{13} = 'cQpv';
    rowHeader{14} = 'cPwi';
    rowHeader{15} = 'cQwi';
    rowHeader{16} = 'cAC';
    rowHeader{17} = 'cClNI';

    Cost = [cV; cPtr; cQtr; cPgen; cQgen; cF; cTr; cSw; cCp; cBatStb; cBatBeta; cPpv; cQpv; cPwi; cQwi; cAC; cClNI];
    sheetName = 'Costos';
    printVarNx1xT(Cost, rowHeader, Header, outFilename, sheetName);
end
