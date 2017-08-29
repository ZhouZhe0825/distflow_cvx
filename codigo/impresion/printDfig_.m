function printDfig_(Header, Var, Data, outFilename)

	if isfield(Var, 'Gen')
		if isfield(Var.Gen, 'Dfig')

			nodDfig = find(Data.Gen.DFIG.I == 1);

			for i = 1:length(nodDfig)
				nod = nodDfig(i);
				Etapas = size(Var.Gen.Dfig.Branch.PIE,2);

				pWi = zeros(1,Etapas);
				qWi = zeros(1,Etapas);
				P = zeros(3,Etapas);
				Q = zeros(3,Etapas);
				l = zeros(3,Etapas);
				r_l = zeros(3,Etapas);
				h_l = zeros(3,Etapas);
				pg = zeros(2,Etapas);
				qg = zeros(2,Etapas);
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
				
				pWi = Var.Gen.Dfig.pWi(nod, :);
				qWi = Var.Gen.Dfig.qWi(nod, :);

				P(1,:) = Var.Gen.Dfig.Branch.PIE(nod, :);
				P(2,:) = Var.Gen.Dfig.Branch.PIF(nod, :);
				P(3,:) = Var.Gen.Dfig.Branch.POR(nod, :);

				Q(1,:) = Var.Gen.Dfig.Branch.QIE(nod, :);
				Q(2,:) = Var.Gen.Dfig.Branch.QIF(nod, :);
				Q(3,:) = Var.Gen.Dfig.Branch.QOR(nod, :);

				l(1,:) = Var.Gen.Dfig.Branch.lIE(nod, :);
				l(2,:) = Var.Gen.Dfig.Branch.lIF(nod, :);
				l(3,:) = Var.Gen.Dfig.Branch.lOR(nod, :);

				r_l(1,:) = Data.Gen.DFIG.rIE(nod, :).*l(1,:);
				r_l(2,:) = Data.Gen.DFIG.rIF(nod, :).*l(2,:);
				r_l(3,:) = Data.Gen.DFIG.rOR(nod, :).*l(3,:);

				v(1,:) =  Var.Gen.Dfig.Bus.vI(nod, :);
				v(2,:) =  Var.Gen.Dfig.Bus.vE(nod, :);
				v(3,:) =  Var.Gen.Dfig.Bus.vF(nod, :);
				v(4,:) =  Var.Gen.Dfig.Bus.vO(nod, :);
				v(5,:) =  Var.Gen.Dfig.Bus.vR(nod, :);
				
				h_l = ((P.^2 + Q.^2) ./ v([1 1 4],:)) ./ (l+eps);
				
				pg(1,:) = Var.Gen.Dfig.Bus.pgE(nod, :);
				pg(2,:) = Var.Gen.Dfig.Bus.pgR(nod, :);

				qg(1,:) = Var.Gen.Dfig.Bus.qgE(nod, :);
				qg(2,:) = Var.Gen.Dfig.Bus.qgR(nod, :);

				pC = Var.Gen.Dfig.Bus.pCF(nod, :);
				qC = Var.Gen.Dfig.Bus.qCF(nod, :);
				
				s(1,:) =  Var.Gen.Dfig.Bus.sF(nod, :);
				s(2,:) =  Var.Gen.Dfig.Bus.sR(nod, :);

				h_s(1,:) = sqrt(pC.^2 + qC.^2)./s(1,:); 
				h_s(2,:) = sqrt(P(3,:).^2 + Q(3,:).^2)./s(2,:);

				cv_s(1,:) = Data.Gen.DFIG.cvF(nod, :).*s(1,:);
				cv_s(2,:) = Data.Gen.DFIG.cvR(nod, :).*s(2,:);
				
				xi(1,:) = Var.Gen.Dfig.Bus.xiF(nod, :);
				xi(2,:) = Var.Gen.Dfig.Bus.xiR(nod, :);

				h_xi(1,:) = (pC.^2 + qC.^2)./xi(1,:); 
				h_xi(2,:) = (P(3,:).^2 + Q(3,:).^2)./xi(2,:); 

				cr_xi(1,:) = Data.Gen.DFIG.crF(nod, :).*xi(1,:);
				cr_xi(2,:) = Data.Gen.DFIG.crR(nod, :).*xi(2,:);

				n_ = Var.Gen.Dfig.Bus.n_Wnd(nod, :);
				P_mec = Var.Gen.Dfig.Bus.P_mecWnd(nod, :);

				
				rowHeader = cell(42,1);		
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
				rowHeader{	21	} = 'qg_5';
				rowHeader{	22	} = 'pC_3';
				rowHeader{	23	} = 'qC_3';
				rowHeader{	24	} = 's_3';
				rowHeader{	25	} = 's_5';
				rowHeader{	26	} = 'Holg s_3';
				rowHeader{	27	} = 'Holg s_5';
				rowHeader{	28	} = 'cv_3 s_3';
				rowHeader{	29	} = 'cv_5 s_5';
				rowHeader{	30	} = 'xi_3';
				rowHeader{	31	} = 'xi_5';
				rowHeader{	32	} = 'Holg xi_3';
				rowHeader{	33	} = 'Holg xi_5';
				rowHeader{	34	} = 'cr_3 xi_3';
				rowHeader{	35	} = 'cr_5 xi_5';
				rowHeader{	36	} = 'v_1';
				rowHeader{	37	} = 'v_2';
				rowHeader{	38	} = 'v_3';
				rowHeader{	39	} = 'v_4';
				rowHeader{	40	} = 'v_5';
				rowHeader{	41	} = 'n_';
				rowHeader{	42	} = 'P_mec';

				Dfig = [pWi;qWi;P;Q;l;h_l;r_l;pg;qg;pC;qC;s;h_s;cv_s;xi;h_xi;cr_xi;v;n_;P_mec];
				sheetName = ['Dfig_' num2str(nod)];
				printVarNx1xT(Dfig, rowHeader, Header, outFilename, sheetName);
			end
		end
	end
end
