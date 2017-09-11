function plotRedPost(Data, Var, t)

	TotalT = matOverTime(Data.Red.Branch.T);
	
	D.Red.Branch.T = TotalT;
	T = repmat(TotalT, [1 1 size(Data.Red.Bus.pCLow,2)]);

	[Var] = VarM2NxN(Var,D);
	n = size(Data.Red.Branch.T,1);

% 	Ncap = find(Data.Red.Bus.Icap == 1);
	Ncons = find(Data.Red.Bus.Icons == 1);
	Nncons = setdiff((1:length(Data.Red.Branch.T)),Ncons);
% 	Ntap = find(Data.Red.Bus.Itap == 1);
% 	Ncarg = find(Data.ClNI.I == 1);
% 	Ndfig = find(Data.Gen.DFIG.I == 1);
% 	Npv = find(Data.Gen.Pv.I == 1);
	Ntras = find(Data.Gen.Tras.I == 1);

	P = Var.Red.Branch.P;
	Q = Var.Red.Branch.Q;
	l = Var.Red.Branch.l;
	v = Var.Red.Bus.v;
	z = Var.Red.Branch.z;
	
	pC = Var.Red.Bus.pC;
	qC = Var.Red.Bus.qC;
	pN = Var.Red.Bus.pN;
	qN = Var.Red.Bus.qN;
	pG = Var.Red.Bus.pG;
	qG = Var.Red.Bus.qG;

	errP = abs(pC - pG) ./ (abs(pN)+eps);
	errQ = abs(qC - qG) ./ (abs(qN)+eps);
	
	PQv = (P.^2 + Q.^2) ./ repmat(v, [1 n 1]).*T;
	lRel = abs(PQv) ./ (abs(l)+eps);
	lRelz = lRel .* z;
	lRelz = (lRelz + permute(lRelz, [2 1 3]));
	
	pMap = [0 .85 .9 .95 .99 2];
	errPD = errP;
	for i = length(pMap):-1:1
		errPD(errP <= pMap(i)) = i*2;
	end
	
% 	Nbat = find(Data.St.Bat.I == 1);

	G = digraph(TotalT);
	EndNodes = G.Edges.EndNodes;
	for i = 1:size(EndNodes,1)
		G.Edges.Weight(i) = lRelz(EndNodes(i,1),EndNodes(i,2),t);
	end
	LWidths = 3*G.Edges.Weight+eps;
	h = plot(G, 'LineStyle', 'none','LineWidth',LWidths);
	layout(h,'layered','Direction','down','Sources',Ntras(1));
	highlight(h, digraph(Var.Red.Branch.z(:,:,t)), 'LineStyle', '-');


	highlight(h,Ncons,'NodeColor','g', 'Marker','^');
	highlight(h,Nncons, 'Marker','o');
	for i = 1:length(T(:,:,t))
		labelnode(h,i,{num2str(i)});
		highlight(h,i,'MarkerSize',errPD(i,:,t));
	end

	title(['T = ', num2str(t)]);
	
	
% 	labelnode(h,Ncap,{'Cap'});
% 	labelnode(h,Ncarg,{'Carga'});
% 	labelnode(h,Ndfig,{'Dfig'});
% 	labelnode(h,Npv,{'Pv'});
% 	labelnode(h,Nbat,{'Bat'});

end