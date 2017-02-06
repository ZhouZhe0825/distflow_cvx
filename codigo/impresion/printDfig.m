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
				r_l = zeros(3,Etapas);
				h_l = zeros(3,Etapas);
				pg = zeros(2,Etapas);
				qg = zeros(3,Etapas);
				pC = zeros(1,Etapas);
				qC = zeros(1,Etapas); 
				s = zeros(2,Etapas);
				h_s = zeros(2,Etapas);
				cv_s = zeros(2,Etapas);
				xi = zeros(2,Etapas);
				h_xi = zeros(2,Etapas);
				cr_xi = zeros(2,Etapas);
				v = zeros(5,Etapas);
				n_ = zeros(1,Etapas);
				P_mec = zeros(1,Etapas); 
				
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

				r_l(1,:) = squeeze(Data.Gen.DFIG.r(1,2,:,i))'.*l(1,:);
				r_l(2,:) = squeeze(Data.Gen.DFIG.r(1,3,:,i))'.*l(2,:);
				r_l(3,:) = squeeze(Data.Gen.DFIG.r(4,5,:,i))'.*l(3,:);

				v(:,:) =  squeeze(Var.Gen.Dfig.Bus.v(:,1,:,i));
				
				h_l = ((P.^2 + Q.^2) ./ v([1 1 4],:)) ./ (l+eps);
				
				pg(:,:) = squeeze(Var.Gen.Dfig.Bus.pg([2 5],1,:,i));
				qg(:,:) = squeeze(Var.Gen.Dfig.Bus.qg([2 3 5],1,:,i));
				pC(:,:) = squeeze(Var.Gen.Dfig.Bus.pC(3,1,:,i));
				qC(:,:) = squeeze(Var.Gen.Dfig.Bus.qC(3,1,:,i));
				
				s(:,:) =  squeeze(Var.Gen.Dfig.Bus.s([3 5],1,:,i));
				h_s(1,:) = sqrt(pC.^2 + qC.^2)./s(1,:); 
				h_s(2,:) = sqrt(P(3,:).^2 + Q(3,:).^2)./s(2,:);
				cv_s(1,:) = squeeze(Data.Gen.DFIG.cv(3,1,:,i))'.*s(1,:);
				cv_s(2,:) = squeeze(Data.Gen.DFIG.cv(5,1,:,i))'.*s(2,:);
				
				xi(:,:) = squeeze(Var.Gen.Dfig.Bus.xi([3 5],1,:,i));
				h_xi(1,:) = (pC.^2 + qC.^2)./xi(1,:); 
				h_xi(2,:) = (P(3,:).^2 + Q(3,:).^2)./xi(2,:); 
				cr_xi(1,:) = squeeze(Data.Gen.DFIG.cr(3,1,:,i))'.*xi(1,:);
				cr_xi(2,:) = squeeze(Data.Gen.DFIG.cr(5,1,:,i))'.*xi(2,:);

				n_ = squeeze(Var.Gen.Dfig.Bus.n_Wnd(1,1,:,i))';
				P_mec = squeeze(Var.Gen.Dfig.Bus.P_mecWnd(1,1,:,i))';

				
				rowHeader = cell(43,1);		
				rowHeader{	1	} = 'pWi';
				rowHeader{	2	} = 'qWi';
				rowHeader{	3	} = 'P_1-2';
				rowHeader{	4	} = 'P_1-3';
				rowHeader{	5	} = 'P_4-5';
				rowHeader{	6	} = 'Q_1-2';
				rowHeader{	7	} = 'Q_1-3';
				rowHeader{	8	} = 'Q_4-5';
				rowHeader{	9	} = 'l_1-2';
				rowHeader{	10	} = 'l_1-3';
				rowHeader{	11	} = 'l_4-5';
				rowHeader{	12	} = 'Holg l_1-2';
				rowHeader{	13	} = 'Holg l_1-3'; 
				rowHeader{	14	} = 'Holg l_4-5'; 
				rowHeader{	15	} = 'r_1-2 l_1-2';
				rowHeader{	16	} = 'r_1-3 l_1-3';
				rowHeader{	17	} = 'r_4-5 l_4-5';
				rowHeader{	18	} = 'pg_2';
				rowHeader{	19	} = 'pg_5';
				rowHeader{	20	} = 'qg_2';
				rowHeader{	21	} = 'qg_3';
				rowHeader{	22	} = 'qg_5';
				rowHeader{	23	} = 'pC_3';
				rowHeader{	24	} = 'qC_3';
				rowHeader{	25	} = 's_3';
				rowHeader{	26	} = 's_5';
				rowHeader{	27	} = 'Holg s_3';
				rowHeader{	28	} = 'Holg s_5';
				rowHeader{	29	} = 'cv_3 s_3';
				rowHeader{	30	} = 'cv_5 s_5';
				rowHeader{	31	} = 'xi_3';
				rowHeader{	32	} = 'xi_5';
				rowHeader{	33	} = 'Holg xi_3';
				rowHeader{	34	} = 'Holg xi_5';
				rowHeader{	35	} = 'cr_3 xi_3';
				rowHeader{	36	} = 'cr_5 xi_5';
				rowHeader{	37	} = 'v_1';
				rowHeader{	38	} = 'v_2';
				rowHeader{	39	} = 'v_3';
				rowHeader{	40	} = 'v_4';
				rowHeader{	41	} = 'v_5';
				rowHeader{	42	} = 'n_';
				rowHeader{	43	} = 'P_mec';

				Dfig = [pWi;qWi;P;Q;l;h_l;r_l;pg;qg;pC;qC;s;h_s;cv_s;xi;h_xi;cr_xi;v;n_;P_mec];
				sheetName = ['Dfig_' num2str(nod)];
				printVarNx1xT(Dfig, rowHeader, Header, outFilename, sheetName);
			end
		end
	end
end
