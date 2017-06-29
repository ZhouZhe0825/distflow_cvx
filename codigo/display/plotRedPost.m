function plotRedPost(Data, Var, t)

% 	Ncap = find(Data.Red.Bus.Icap == 1);
% 	Ncons = find(Data.Red.Bus.Icons == 1);
% 	Nncons = setdiff((1:length(Data.Red.Branch.T)),Ncons);
% 	Ntap = find(Data.Red.Bus.Itap == 1);
% 	Ncarg = find(Data.ClNI.I == 1);
% 	Ndfig = find(Data.Gen.DFIG.I == 1);
% 	Npv = find(Data.Gen.Pv.I == 1);
	Ntras = find(Data.Gen.Tras.I == 1);
% 	Nbat = find(Data.St.Bat.I == 1);



	G = digraph(Data.Red.Branch.T(:,:,t));
	h = plot(G, 'LineStyle', 'none');
	layout(h,'layered','Direction','down','Sources',Ntras(1));
	highlight(h, digraph(Var.Red.Branch.z(:,:,t)), 'LineStyle', '-');
% 	highlight(h,Ncons,'NodeColor','g', 'Marker','^','MarkerSize',4);
% 	highlight(h,Nncons, 'Marker','none');
    for i = 1:length(Data.Red.Branch.T(:,:,t))
        labelnode(h,i,{num2str(i)});
    end
% 	labelnode(h,Ncap,{'Cap'});
% 	labelnode(h,Ncarg,{'Carga'});
% 	labelnode(h,Ndfig,{'Dfig'});
% 	labelnode(h,Npv,{'Pv'});
% 	labelnode(h,Nbat,{'Bat'});

end