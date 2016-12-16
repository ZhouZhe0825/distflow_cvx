function printDfig(Header, Var, Data, outFilename)

	if isfield(Var, 'Gen')
		if isfield(Var.Gen, 'Dfig')

			nodDfig = find(matOverTime(Data.Gen.DFIG.I) == 1);

			for i = 1:length(nodDfig)
				nod = nodDfig(i);
				Etapas = size(Var.Gen.Dfig.Branch.P,3);

				pWi = zeros(1,Etapas);
				qWi = zeros(1,Etapas);
				P = zeros(3,Etapas);
				Q = zeros(3,Etapas);
				l = zeros(3,Etapas);
				pg = zeros(2,Etapas);
				qg = zeros(3,Etapas);
				pC = zeros(1,Etapas);
				s = zeros(2,Etapas);
				xi = zeros(2,Etapas);
				v = zeros(5,Etapas);
				
				pWi(:,:) = squeeze(Var.Gen.Dfig.pWi(nod,1,:));
				qWi(:,:) = squeeze(Var.Gen.Dfig.qWi(nod,1,:));

				P(1,:) = squeeze(Var.Gen.Dfig.Branch.P(1,2,:,i));
				P(2,:) = squeeze(Var.Gen.Dfig.Branch.P(1,3,:,i));
				P(3,:) = squeeze(Var.Gen.Dfig.Branch.P(4,5,:,i));

				Q(1,:) = squeeze(Var.Gen.Dfig.Branch.Q(1,2,:,i));
				Q(2,:) = squeeze(Var.Gen.Dfig.Branch.Q(1,3,:,i));
				Q(3,:) = squeeze(Var.Gen.Dfig.Branch.Q(4,5,:,i));

				l(1,:) = squeeze(Var.Gen.Dfig.Branch.l(1,2,:,i));
				l(2,:) = squeeze(Var.Gen.Dfig.Branch.l(1,3,:,i));
				l(3,:) = squeeze(Var.Gen.Dfig.Branch.l(4,5,:,i));

				pg(:,:) = squeeze(Var.Gen.Dfig.Bus.pg([2 5],1,:,i));
				qg(:,:) = squeeze(Var.Gen.Dfig.Bus.qg([2 3 5],1,:,i));
				pC(:,:) = squeeze(Var.Gen.Dfig.Bus.pC(3,1,:,i));
				s(:,:) =  squeeze(Var.Gen.Dfig.Bus.s([3 5],1,:,i));
				xi(:,:) = squeeze(Var.Gen.Dfig.Bus.xi([3 5],1,:,i));
				v(:,:) =  squeeze(Var.Gen.Dfig.Bus.v(:,1,:,i));

				rowHeader = cell(26,1);
				rowHeader{1} = 'pWi';
				rowHeader{2} = 'qWi';
				rowHeader{3} = 'P_1-2';
				rowHeader{4} = 'P_1-3';
				rowHeader{5} = 'P_4-5';
				rowHeader{6} = 'Q_1-2';
				rowHeader{7} = 'Q_1-3';
				rowHeader{8} = 'Q_4-5';
				rowHeader{9} = 'l_1-2';
				rowHeader{10} = 'l_1-3';
				rowHeader{11} = 'l_4-5';
				rowHeader{12} = 'pg_2';
				rowHeader{13} = 'pg_5';
				rowHeader{14} = 'qg_2';
				rowHeader{15} = 'qg_3';
				rowHeader{16} = 'qg_5';
				rowHeader{17} = 'pC_3';
				rowHeader{18} = 's_3';
				rowHeader{19} = 's_5';
				rowHeader{20} = 'xi_3';
				rowHeader{21} = 'xi_5';
				rowHeader{22} = 'v_1';
				rowHeader{23} = 'v_2';
				rowHeader{24} = 'v_3';
				rowHeader{25} = 'v_4';
				rowHeader{26} = 'v_5';

				Dfig = [pWi;qWi;P;Q;l;pg;qg;pC;s;xi;v];
				sheetName = ['Dfig_' num2str(nod)];
				printVarNx1xT(Dfig, rowHeader, Header, outFilename, sheetName);
			end
		end
	end
end
