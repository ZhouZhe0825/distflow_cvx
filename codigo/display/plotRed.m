function plotRed(Data)

	Ncap = find(Data.Red.Bus.Icap == 1);
	Ncons = find(Data.Red.Bus.Icons == 1);
	Nncons = setdiff((1:length(Data.Red.Branch.T)),Ncons);
	[NtapI, NtapJ, ~] = find(triu(Data.Red.Branch.Itap) == 1);
	Ncarg = find(Data.ClNI.I == 1);
	Ndfig = find(Data.Gen.DFIG.I == 1);
	Npv = find(Data.Gen.Pv.I == 1);
	Ntras = find(Data.Gen.Tras.I == 1);
	Nbat = find(Data.St.Bat.I == 1);

    


	G = graph(Data.Red.Branch.T);
	h = plot(G);
	labelnode(h,Ntras,{'Tras'});
	layout(h,'layered','Direction','down','Sources',Ntras(1));
	highlight(h, graph(Data.Red.Branch.Tswitches), 'EdgeColor', 'r', 'LineStyle', '--', 'Linewidth',1.5);
	highlight(h,Ncons,'NodeColor','g', 'Marker','^','MarkerSize',4);
    for i = 1:length(Data.Red.Branch.T(:,:))
        labelnode(h,i,{num2str(i)});
    end
	highlight(h,Nncons, 'Marker','none');
	highlight(h, graph(Data.Red.Branch.Tswitches), 'EdgeColor', 'r', 'LineStyle', '--', 'Linewidth',1.5);
    labeledge(h,NtapI,NtapJ,repmat({'Tap'},[length(NtapI) 1]))    
	labelnode(h,Ncap,{'Cap'});
	labelnode(h,Ncarg,{'Carga'});
	labelnode(h,Ndfig,{'Dfig'});
	labelnode(h,Npv,{'Pv'});
	labelnode(h,Nbat,{'Bat'});

end