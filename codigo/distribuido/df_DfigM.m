function [Var, opt, status] = df_Dfig(Data, Config, DistrInfo)

%% Inicializacion
VertI = VertIMat(Data.Red.Branch.T);
VertJ = VertJMat(Data.Red.Branch.T);
OutBr = VertI';
InBr = VertJ';

n = size(VertI,2);
m = size(VertI,1);

G = find(Data.Gen.Tras.I == 1);
NnoG = setdiff((1:n),G);

Data.Red.Branch.Tup = triu(Data.Red.Branch.T);
[rowTup, colTup, ~] = find(Data.Red.Branch.Tup == 1);
Tind = find(Data.Red.Branch.T == 1);
Tupmind = sub2ind(size(Data.Red.Branch.T), rowTup, colTup);
Tdownmind = sub2ind(size(Data.Red.Branch.T), colTup, rowTup);
Tupm = zeros(size(Tupmind));
Tdownm = zeros(size(Tdownmind));
for i=1:length(Tupmind)
	Tupm(i) = find(Tind == Tupmind(i));
	Tdownm(i) = find(Tind == Tdownmind(i));
end

tnnLow = (1 + Data.Red.Bus.TapLow.*Data.Red.Bus.Ntr);
tnnTop = (1 + Data.Red.Bus.TapTop.*Data.Red.Bus.Ntr);

NcpCapL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow;
NcpCapT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop;

NcpvL = Data.Red.Bus.Ncp.*(Data.Red.Bus.uLow).^2;
NcpvT = Data.Red.Bus.Ncp.*(Data.Red.Bus.uTop).^2;

NcpCapLvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uLow).^2;
NcpCapTvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uTop).^2;
NcpCapLvT = Data.Red.Bus.Ncp.*Data.Red.Bus.CapLow.*(Data.Red.Bus.uTop).^2;
NcpCapTvL = Data.Red.Bus.Ncp.*Data.Red.Bus.CapTop.*(Data.Red.Bus.uLow).^2;

% Eolico
indWn = DistrInfo.Bus;
NindWn = setdiff((1:n),indWn);
lenWN = length(indWn);

P_mecSigWnd = abs(sign(Data.Gen.DFIG.P_mec));
P_mecWnd = Data.Gen.DFIG.P_mec;
n_Wnd = Data.Gen.DFIG.n_;

%% Modelo programacion matematica

tic
cvx_begin quiet

	for i = 1: size(Config.SubP,1)
		cvx_solver_settings(Config.SubP{i,1}, Config.SubP{i,2});
	end

	cvx_precision high
	%% Declaracion de variables y expresiones
	variable pWi(lenWN,Config.Etapas);
	variable qWi(lenWN,Config.Etapas);

	expression tfopt_expr(1,Config.Etapas); 
	expression tfopt_virt(1,Config.Etapas); 
	expression tfopt_conv(1,Config.Etapas); 
	expression fopt_expr; 

	% Variables de generadores eolico
	variable cqWi(lenWN, Config.Etapas);

	variable PdfigIE(lenWN,Config.Etapas);
	variable PdfigIF(lenWN,Config.Etapas);
	variable PdfigOR(lenWN,Config.Etapas);

	variable QdfigIE(lenWN,Config.Etapas);
	variable QdfigIF(lenWN,Config.Etapas);
	variable QdfigOR(lenWN,Config.Etapas);

	variable ldfigIE(lenWN,Config.Etapas);
	variable ldfigIF(lenWN,Config.Etapas);
	variable ldfigOR(lenWN,Config.Etapas);

	variable vdfigI(lenWN,Config.Etapas);
	variable vdfigE(lenWN,Config.Etapas);
	variable vdfigF(lenWN,Config.Etapas);
	variable vdfigO(lenWN,Config.Etapas);
	variable vdfigR(lenWN,Config.Etapas);

	variable pWigdfigE(lenWN,Config.Etapas);
	variable pWigdfigR(lenWN,Config.Etapas);

	variable qWigdfigE(lenWN,Config.Etapas);
	variable qWigdfigR(lenWN,Config.Etapas);

	variable pCdfigF(lenWN,Config.Etapas);
	variable qCdfigF(lenWN,Config.Etapas);

	variable sdfigF(lenWN,Config.Etapas);
	variable sdfigR(lenWN,Config.Etapas);

	variable xidfigF(lenWN,Config.Etapas);
	variable xidfigR(lenWN,Config.Etapas);

	expression lQoLdfigIE(lenWN,Config.Etapas,3);
	expression lQoLdfigIF(lenWN,Config.Etapas,3);
	expression lQoLdfigOR(lenWN,Config.Etapas,3);

	expression sNormdfigF(lenWN,Config.Etapas,2);
	expression sNormdfigR(lenWN,Config.Etapas,2);

	expression PQNormdfigIE(lenWN,Config.Etapas,2);
	expression PQNormdfigIF(lenWN,Config.Etapas,2);


	tfopt_expr = ...
		sum(Data.Cost.rhopWi(indWn) .* pWi) ...
		+ sum(cqWi) ...
	;

	tfopt_conv = 1/(2*DistrInfo.Gama) * ...
				sum(norms(pWi - DistrInfo.Dfig.pWi(indWn,:),2,2) ...
				+ norms(qWi - DistrInfo.Dfig.qWi(indWn,:),2,2),1);

	tfopt_virt = sum(DistrInfo.muT(indWn,:) .* pWi + DistrInfo.lambdaT(indWn,:) .* qWi,1);

	fopt_expr = sum(tfopt_expr) + sum(tfopt_virt) + sum(tfopt_conv);

	cqWi >= - Data.Cost.rhomqWi(indWn) .* qWi;
	cqWi >= Data.Cost.rhoMqWi(indWn) .* qWi;

	pWi == - (PdfigIE + PdfigIF);
	qWi == - (QdfigIE + QdfigIF);

	% Modelo de Red interna
	vdfigI == DistrInfo.Dfig.v(indWn,1,:);

	PdfigIE == Data.Gen.DFIG.rIE .* ldfigIE - pWigdfigE;
	PdfigIF == Data.Gen.DFIG.rIF .* ldfigIF + pCdfigF;
	PdfigOR == Data.Gen.DFIG.rOR .* ldfigOR - pWigdfigR;

	QdfigIE == Data.Gen.DFIG.xIE .* ldfigIE - qWigdfigE;
	QdfigIF == Data.Gen.DFIG.xIF .* ldfigIF - qCdfigF;
	QdfigOR == n_Wnd .* Data.Gen.DFIG.xOR .* ldfigOR - qWigdfigR;

	vdfigE == vdfigI - 2*(Data.Gen.DFIG.rIE .* PdfigIE + Data.Gen.DFIG.xIE .* QdfigIE) + (Data.Gen.DFIG.rIE.^2 + Data.Gen.DFIG.xIE.^2) .* ldfigIE;
	vdfigF == vdfigI - 2*(Data.Gen.DFIG.rIF .* PdfigIF + Data.Gen.DFIG.xIF .* QdfigIF) + (Data.Gen.DFIG.rIF.^2 + Data.Gen.DFIG.xIF.^2) .* ldfigIF;
	vdfigR == vdfigO - 2*(Data.Gen.DFIG.rOR .* PdfigOR + n_Wnd.*Data.Gen.DFIG.xOR .* QdfigOR) + (Data.Gen.DFIG.rOR.^2 + n_Wnd.^2 .* Data.Gen.DFIG.xOR.^2) .* ldfigOR;

	% Corriente
	lQoLdfigIE(:,:,1) = 2*PdfigIE;
	lQoLdfigIE(:,:,2) = 2*QdfigIE;
	lQoLdfigIE(:,:,3) = ldfigIE - vdfigI;
	norms(lQoLdfigIE,2,3) - (ldfigIE + vdfigI) <= 0;

	lQoLdfigIF(:,:,1) = 2*PdfigIF;
	lQoLdfigIF(:,:,2) = 2*QdfigIF;
	lQoLdfigIF(:,:,3) = ldfigIF - vdfigI;
	norms(lQoLdfigIF,2,3) - (ldfigIF + vdfigI) <= 0;

	lQoLdfigOR(:,:,1) = 2*PdfigOR;
	lQoLdfigOR(:,:,2) = 2*QdfigOR;
	lQoLdfigOR(:,:,3) = ldfigOR - vdfigO;
	norms(lQoLdfigOR,2,3) - (ldfigOR + vdfigO) <= 0;

	Data.Gen.DFIG.lTopIF >= ldfigIF;
	Data.Gen.DFIG.lTopOR >= ldfigOR;

	(Data.Gen.DFIG.lTopIE - ldfigIE).*P_mecSigWnd >= 0;
	(Data.Gen.DFIG.sTopF - sdfigF).*P_mecSigWnd >= 0;
	(Data.Gen.DFIG.sTopR - sdfigR).*P_mecSigWnd >= 0;
	(Data.Gen.DFIG.xiTopF - xidfigF).*P_mecSigWnd >= 0;
	(Data.Gen.DFIG.xiTopR - xidfigR).*P_mecSigWnd >= 0;

	vdfigE >= Data.Gen.DFIG.uLowE.^2;
	vdfigE <= Data.Gen.DFIG.uTopE.^2;

	vdfigF >= Data.Gen.DFIG.uLowF.^2;
	vdfigF <= Data.Gen.DFIG.uTopF.^2;

	PQNormdfigIE(:,:,1) = PdfigIE;
	PQNormdfigIE(:,:,2) = QdfigIE;
	Data.Gen.DFIG.PQnormIE >= norms(PQNormdfigIE,2,3);

	PQNormdfigIF(:,:,1) = PdfigIF;
	PQNormdfigIF(:,:,2) = QdfigIF;
	Data.Gen.DFIG.PQnormIF >= norms(PQNormdfigIF,2,3);

	sNormdfigF(:,:,1) = pCdfigF;
	sNormdfigF(:,:,2) = qCdfigF;
	sdfigF >= norms(sNormdfigF,2,3);
	xidfigF >= pCdfigF.^2 + qCdfigF.^2;

	sNormdfigR(:,:,1) = PdfigOR;
	sNormdfigR(:,:,2) = QdfigOR;
	sdfigR >= norms(sNormdfigR,2,3);
	xidfigR >= PdfigOR.^2 + QdfigOR.^2;

	pCdfigF == PdfigOR ...
		+ (Data.Gen.DFIG.cvR .* sdfigR + Data.Gen.DFIG.crR .* xidfigR) ...
		+ (Data.Gen.DFIG.cvF .* sdfigF + Data.Gen.DFIG.crF .* xidfigF);
	vdfigR == n_Wnd.^2*(Data.Gen.DFIG.N_er^2).*vdfigE;
	pWigdfigE == P_mecWnd ./ (1-n_Wnd); %TODO es igual
	pWigdfigR == -n_Wnd.*pWigdfigE;
	qWigdfigR == n_Wnd.*(qWigdfigE);


	fopt_expr = sum(tfopt_expr + tfopt_virt + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT
Var.Gen.Dfig.pWi = zeros(n,Config.Etapas);
Var.Gen.Dfig.qWi = zeros(n,Config.Etapas);
Var.Gen.Dfig.pWi(indWn,:) = pWi;
Var.Gen.Dfig.qWi(indWn,:) = qWi;

Var.Gen.Dfig.Branch.PIE = PdfigIE';
Var.Gen.Dfig.Branch.PIF = PdfigIF';
Var.Gen.Dfig.Branch.POR = PdfigOR';

Var.Gen.Dfig.Branch.QIE = QdfigIE';
Var.Gen.Dfig.Branch.QIF = QdfigIF';
Var.Gen.Dfig.Branch.QOR = QdfigOR';

Var.Gen.Dfig.Branch.lIE = ldfigIE';
Var.Gen.Dfig.Branch.lIF = ldfigIF';
Var.Gen.Dfig.Branch.lOR = ldfigOR';

Var.Gen.Dfig.Bus.vI = vdfigI';
Var.Gen.Dfig.Bus.vE = vdfigE';
Var.Gen.Dfig.Bus.vF = vdfigF';
Var.Gen.Dfig.Bus.vO = vdfigO';
Var.Gen.Dfig.Bus.vR = vdfigR';

Var.Gen.Dfig.Bus.pCF = pCdfigF';

Var.Gen.Dfig.Bus.qCF = qCdfigF';

Var.Gen.Dfig.Bus.pgE = pWigdfigE';
Var.Gen.Dfig.Bus.pgR = pWigdfigR';

Var.Gen.Dfig.Bus.qgE = qWigdfigE';
Var.Gen.Dfig.Bus.qgR = qWigdfigR';

Var.Gen.Dfig.Bus.sF = sdfigF';
Var.Gen.Dfig.Bus.sR = sdfigR';

Var.Gen.Dfig.Bus.xiF = xidfigF';
Var.Gen.Dfig.Bus.xiR = xidfigR';

Var.Gen.Dfig.Bus.n_Wnd = n_Wnd;
Var.Gen.Dfig.Bus.P_mecWnd = P_mecWnd;

Var.Red.Bus.pG = Var.Gen.Dfig.pWi;
Var.Red.Bus.qG = Var.Gen.Dfig.qWi;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_virt);
opt(1,3) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end