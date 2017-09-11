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
	expression tfopt_mu(Config.Etapas,1); 
	expression tfopt_lambda(Config.Etapas,1); 
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
		sum(Data.Cost.rhopWi(indWn) .* pWi,1) ...
		+ sum(cqWi,1) ...
	;

	tfopt_conv = 1/(2*DistrInfo.Gama) * ...
				(norms(pWi - DistrInfo.Dfig.pWi(indWn,:),2,1) ...
				+ norms(qWi - DistrInfo.Dfig.qWi(indWn,:),2,2));

	tfopt_mu = sum(DistrInfo.muT(indWn,:) .* pWi,1);
	tfopt_lambda = sum(DistrInfo.lambdaT(indWn,:) .* qWi,1);

	cqWi >= - Data.Cost.rhomqWi(indWn) .* qWi;
	cqWi >= Data.Cost.rhoMqWi(indWn) .* qWi;

	pWi == - (PdfigIE + PdfigIF);
	qWi == - (QdfigIE + QdfigIF);

	% Modelo de Red interna
	vdfigI == DistrInfo.Dfig.v(indWn,1,:);

	PdfigIE == Data.Gen.DFIG.rIE(indWn,:) .* ldfigIE - pWigdfigE;
	PdfigIF == Data.Gen.DFIG.rIF(indWn,:) .* ldfigIF + pCdfigF;
	PdfigOR == Data.Gen.DFIG.rOR(indWn,:) .* ldfigOR - pWigdfigR;

	QdfigIE == Data.Gen.DFIG.xIE(indWn,:) .* ldfigIE - qWigdfigE;
	QdfigIF == Data.Gen.DFIG.xIF(indWn,:) .* ldfigIF - qCdfigF;
	QdfigOR == n_Wnd(indWn,:) .* Data.Gen.DFIG.xOR(indWn,:) .* ldfigOR - qWigdfigR;

	vdfigE == vdfigI - 2*(Data.Gen.DFIG.rIE(indWn,:) .* PdfigIE + Data.Gen.DFIG.xIE(indWn,:) .* QdfigIE) + (Data.Gen.DFIG.rIE(indWn,:).^2 + Data.Gen.DFIG.xIE(indWn,:).^2) .* ldfigIE;
	vdfigF == vdfigI - 2*(Data.Gen.DFIG.rIF(indWn,:) .* PdfigIF + Data.Gen.DFIG.xIF(indWn,:) .* QdfigIF) + (Data.Gen.DFIG.rIF(indWn,:).^2 + Data.Gen.DFIG.xIF(indWn,:).^2) .* ldfigIF;
	vdfigR == vdfigO - 2*(Data.Gen.DFIG.rOR(indWn,:) .* PdfigOR + n_Wnd(indWn,:).*Data.Gen.DFIG.xOR(indWn,:) .* QdfigOR) + (Data.Gen.DFIG.rOR(indWn,:).^2 + n_Wnd(indWn,:).^2 .* Data.Gen.DFIG.xOR(indWn,:).^2) .* ldfigOR;

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

	Data.Gen.DFIG.lTopIF(indWn,:) >= ldfigIF;
	Data.Gen.DFIG.lTopOR(indWn,:) >= ldfigOR;

	(Data.Gen.DFIG.lTopIE(indWn,:)).*P_mecSigWnd(indWn,:) - ldfigIE >= 0;
	(Data.Gen.DFIG.sTopF(indWn,:)).* P_mecSigWnd(indWn,:) - sdfigF >= 0;
	(Data.Gen.DFIG.sTopR(indWn,:)) .* P_mecSigWnd(indWn,:) - sdfigR >= 0;
	(Data.Gen.DFIG.xiTopF(indWn,:)).*P_mecSigWnd(indWn,:) - xidfigF >= 0;
	(Data.Gen.DFIG.xiTopR(indWn,:)).*P_mecSigWnd(indWn,:) - xidfigR >= 0;

	vdfigE >= Data.Gen.DFIG.uLowE(indWn,:).^2;
	vdfigE <= Data.Gen.DFIG.uTopE(indWn,:).^2;

	vdfigF >= Data.Gen.DFIG.uLowF(indWn,:).^2;
	vdfigF <= Data.Gen.DFIG.uTopF(indWn,:).^2;

	PQNormdfigIE(:,:,1) = PdfigIE;
	PQNormdfigIE(:,:,2) = QdfigIE;
	Data.Gen.DFIG.PQnormIE(indWn,:) >= norms(PQNormdfigIE,2,3);

	PQNormdfigIF(:,:,1) = PdfigIF;
	PQNormdfigIF(:,:,2) = QdfigIF;
	Data.Gen.DFIG.PQnormIF(indWn,:) >= norms(PQNormdfigIF,2,3);

	sNormdfigF(:,:,1) = pCdfigF;
	sNormdfigF(:,:,2) = qCdfigF;
	sdfigF >= norms(sNormdfigF,2,3);
	xidfigF >= pCdfigF.^2 + qCdfigF.^2;

	sNormdfigR(:,:,1) = PdfigOR;
	sNormdfigR(:,:,2) = QdfigOR;
	sdfigR >= norms(sNormdfigR,2,3);
	xidfigR >= PdfigOR.^2 + QdfigOR.^2;

	pCdfigF == PdfigOR ...
		+ (Data.Gen.DFIG.cvR(indWn,:) .* sdfigR + Data.Gen.DFIG.crR(indWn,:) .* xidfigR) ...
		+ (Data.Gen.DFIG.cvF(indWn,:) .* sdfigF + Data.Gen.DFIG.crF(indWn,:) .* xidfigF);
	vdfigR == n_Wnd(indWn,:).^2.*(repmat(Data.Gen.DFIG.N_er(indWn,:).^2,[1 Config.Etapas])).*vdfigE;
	pWigdfigE == P_mecWnd(indWn,:) ./ (1-n_Wnd(indWn,:)); %TODO es igual
	pWigdfigR == -n_Wnd(indWn,:).*pWigdfigE;
	qWigdfigR == n_Wnd(indWn,:).*(qWigdfigE);


	fopt_expr = sum(tfopt_expr + tfopt_mu + tfopt_lambda + tfopt_conv);
	minimize fopt_expr

cvx_end
toc

%% Construccion de la estructura de solucion
% pasaje a NxNxT
Var.Gen.Dfig.pWi = zeros(n,Config.Etapas);
Var.Gen.Dfig.qWi = zeros(n,Config.Etapas);
Var.Gen.Dfig.cqWi = zeros(n,Config.Etapas);
Var.Gen.Dfig.pWi(indWn,:) = pWi;
Var.Gen.Dfig.qWi(indWn,:) = qWi;
Var.Gen.Dfig.cqWi(indWn,:) = cqWi;

Var.Gen.Dfig.Branch.PIE = zeros(n,Config.Etapas);
Var.Gen.Dfig.Branch.PIE(indWn,:) = PdfigIE;
Var.Gen.Dfig.Branch.PIF = zeros(n,Config.Etapas);
Var.Gen.Dfig.Branch.PIF(indWn,:) = PdfigIF;
Var.Gen.Dfig.Branch.POR = zeros(n,Config.Etapas);
Var.Gen.Dfig.Branch.POR(indWn,:) = PdfigOR;

Var.Gen.Dfig.Branch.QIE = zeros(n,Config.Etapas);
Var.Gen.Dfig.Branch.QIE(indWn,:) = QdfigIE;
Var.Gen.Dfig.Branch.QIF = zeros(n,Config.Etapas);
Var.Gen.Dfig.Branch.QIF(indWn,:) = QdfigIF;
Var.Gen.Dfig.Branch.QOR = zeros(n,Config.Etapas);
Var.Gen.Dfig.Branch.QOR(indWn,:) = QdfigOR;

Var.Gen.Dfig.Branch.lIE = zeros(n,Config.Etapas);
Var.Gen.Dfig.Branch.lIE(indWn,:) = ldfigIE;
Var.Gen.Dfig.Branch.lIF = zeros(n,Config.Etapas);
Var.Gen.Dfig.Branch.lIF(indWn,:) = ldfigIF;
Var.Gen.Dfig.Branch.lOR = zeros(n,Config.Etapas);
Var.Gen.Dfig.Branch.lOR(indWn,:) = ldfigOR;

Var.Gen.Dfig.Bus.vI = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.vI(indWn,:) = vdfigI;
Var.Gen.Dfig.Bus.vE = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.vE(indWn,:) = vdfigE;
Var.Gen.Dfig.Bus.vF = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.vF(indWn,:) = vdfigF;
Var.Gen.Dfig.Bus.vO = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.vO(indWn,:) = vdfigO;
Var.Gen.Dfig.Bus.vR = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.vR(indWn,:) = vdfigR;

Var.Gen.Dfig.Bus.pCF = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.pCF(indWn,:) = pCdfigF;

Var.Gen.Dfig.Bus.qCF = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.qCF(indWn,:) = qCdfigF;

Var.Gen.Dfig.Bus.pgE = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.pgE(indWn,:) = pWigdfigE;
Var.Gen.Dfig.Bus.pgR = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.pgR(indWn,:) = pWigdfigR;

Var.Gen.Dfig.Bus.qgE = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.qgE(indWn,:) = qWigdfigE;
Var.Gen.Dfig.Bus.qgR = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.qgR(indWn,:) = qWigdfigR;

Var.Gen.Dfig.Bus.sF = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.sF(indWn,:) = sdfigF;
Var.Gen.Dfig.Bus.sR = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.sR(indWn,:) = sdfigR;

Var.Gen.Dfig.Bus.xiF = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.xiF(indWn,:) = xidfigF;
Var.Gen.Dfig.Bus.xiR = zeros(n,Config.Etapas);
Var.Gen.Dfig.Bus.xiR(indWn,:) = xidfigR;

Var.Gen.Dfig.Bus.n_Wnd = n_Wnd;
Var.Gen.Dfig.Bus.P_mecWnd = P_mecWnd;

Var.Red.Bus.pG = Var.Gen.Dfig.pWi;
Var.Red.Bus.qG = Var.Gen.Dfig.qWi;

opt(1,1) = sum(tfopt_expr);
opt(1,2) = sum(tfopt_mu);
opt(1,3) = sum(tfopt_lambda);
opt(1,4) = sum(tfopt_conv);

status = cvx_status;
cvx_clear
end