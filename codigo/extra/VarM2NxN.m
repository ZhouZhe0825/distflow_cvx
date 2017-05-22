function [VarNxN] = VarM2NxN(VarM, Data)

	VertI = VertIMat(Data.Red.Branch.T);
	VertJ = VertJMat(Data.Red.Branch.T);
	n = size(Data.Red.Branch.T,1);
	Etapas = size(VarM.Red.Branch.P,2);

	VarNxN.Red.Branch.P	 = MxT2NxNxT(VertI,VertJ,	VarM.Red.Branch.P	);
	VarNxN.Red.Branch.Q	 = MxT2NxNxT(VertI,VertJ,	VarM.Red.Branch.Q	);
	VarNxN.Red.Branch.l	 = MxT2NxNxT(VertI,VertJ,	VarM.Red.Branch.l	);
	VarNxN.Red.Branch.z	 = MxT2NxNxT(VertI,VertJ,	VarM.Red.Branch.z	);
	VarNxN.Red.Branch.y	 = MxT2NxNxT(VertI,VertJ,	VarM.Red.Branch.y	);
	VarNxN.Red.Bus.w	 = permute(full(	VarM.Red.Bus.w	), [1 3 2]);

	VarNxN.Red.Bus.v	 = permute(full(	VarM.Red.Bus.v	), [1 3 2]);
	VarNxN.Red.Bus.cDv	 = permute(full(	VarM.Red.Bus.cDv	), [1 3 2]);
	VarNxN.Red.Bus.nn	 = permute(full(	VarM.Red.Bus.nn	), [1 3 2]);
	VarNxN.Red.Bus.nv	 = permute(full(	VarM.Red.Bus.nv	), [1 3 2]);
	VarNxN.Red.Bus.Ntr	 = permute(full(	VarM.Red.Bus.Ntr	), [1 3 2]);

	VarNxN.Red.Bus.pC	 = permute(full(	VarM.Red.Bus.pC	), [1 3 2]);
	VarNxN.Red.Bus.qC	 = permute(full(	VarM.Red.Bus.qC	), [1 3 2]);

	VarNxN.Red.Bus.pN	 = permute(full(	VarM.Red.Bus.pN	), [1 3 2]);
	VarNxN.Red.Bus.qN	 = permute(full(	VarM.Red.Bus.qN	), [1 3 2]);

	VarNxN.Red.Bus.pG	 = permute(full(	VarM.Red.Bus.pG	), [1 3 2]);
	VarNxN.Red.Bus.qG	 = permute(full(	VarM.Red.Bus.qG	), [1 3 2]);

	VarNxN.Red.Bus.qCp	 = permute(full(	VarM.Red.Bus.qCp	), [1 3 2]);
	VarNxN.Red.Bus.Ncp	 = permute(full(	VarM.Red.Bus.Ncp	), [1 3 2]);

	VarNxN.Red.Bus.PTras	 = permute(full(	VarM.Red.Bus.PTras	), [1 3 2]);
	VarNxN.Red.Bus.QTras	 = permute(full(	VarM.Red.Bus.QTras	), [1 3 2]);

	VarNxN.ClRes.pCApp	 = permute(full(	VarM.ClRes.pCApp	), [1 4 2 3]);
	VarNxN.ClRes.qCApp	 = permute(full(	VarM.ClRes.qCApp	), [1 4 2 3]);
	VarNxN.ClRes.pC	 = permute(full(	VarM.ClRes.pC	), [1 3 2]);
	VarNxN.ClRes.qC	 = permute(full(	VarM.ClRes.qC	), [1 3 2]);

	if isfield(VarM, 'Dual')
		VarNxN.Dual.dPn	 = 	VarM.Dual.dPn;
		VarNxN.Dual.dQn	 = 	VarM.Dual.dQn;
	end

	% Aire Acondicionado
	if isfield(VarM, 'ClRes')
		if isfield(VarM.ClRes, 'Tvar')
			VarNxN.ClRes.Tvar	 = permute(full(	VarM.ClRes.Tvar	), [1 3 2]);
		end
	end

	% Baterias
	if isfield(VarM, 'St')
		if isfield(VarM.St, 'Bat')

			VarNxN.St.Bat.pStb	 = permute(full(	VarM.St.Bat.pStb	), [1 3 2]);
			VarNxN.St.Bat.pStgb	 = permute(full(	VarM.St.Bat.pStgb	), [1 3 2]);
			VarNxN.St.Bat.qStb	 = permute(full(	VarM.St.Bat.qStb	), [1 3 2]);
			VarNxN.St.Bat.sStb	 = permute(full(	VarM.St.Bat.sStb	), [1 3 2]);
			VarNxN.St.Bat.xiStb	 = permute(full(	VarM.St.Bat.xiStb	), [1 3 2]);
			VarNxN.St.Bat.EStb	 = permute(full(	VarM.St.Bat.EStb	), [1 3 2]);
		end
	end

	% Cargas No interrumpibles
	if isfield(VarM, 'ClNI')

		VarNxN.ClNI.pC	 = permute(full(	VarM.ClNI.pC	), [1 3 2]);
		VarNxN.ClNI.qC	 = permute(full(	VarM.ClNI.qC	), [1 3 2]);
		VarNxN.ClNI.on	 = permute(full(	VarM.ClNI.on	), [1 3 2]);
		VarNxN.ClNI.start	 = permute(full(	VarM.ClNI.start	), [1 3 2]);

	end

	% Eolico
	if isfield(VarM, 'Gen')
		if isfield(VarM.Gen, 'Dfig')
			lenWN = size(VarM.Gen.Dfig.Branch.PIE,2);

			VarNxN.Gen.Dfig.pWi	 = permute(full(	VarM.Gen.Dfig.pWi	), [1 3 2]);
			VarNxN.Gen.Dfig.qWi	 = permute(full(	VarM.Gen.Dfig.qWi	), [1 3 2]);

			VarNxN.Gen.Dfig.Branch.PIE	 = permute(full(	VarM.Gen.Dfig.Branch.PIE	), [1 3 2]);
			VarNxN.Gen.Dfig.Branch.PIF	 = permute(full(	VarM.Gen.Dfig.Branch.PIF	), [1 3 2]);
			VarNxN.Gen.Dfig.Branch.POR	 = permute(full(	VarM.Gen.Dfig.Branch.POR	), [1 3 2]);

			VarNxN.Gen.Dfig.Branch.QIE	 = permute(full(	VarM.Gen.Dfig.Branch.QIE	), [1 3 2]);
			VarNxN.Gen.Dfig.Branch.QIF	 = permute(full(	VarM.Gen.Dfig.Branch.QIF	), [1 3 2]);
			VarNxN.Gen.Dfig.Branch.QOR	 = permute(full(	VarM.Gen.Dfig.Branch.QOR	), [1 3 2]);

			VarNxN.Gen.Dfig.Branch.lIE	 = permute(full(	VarM.Gen.Dfig.Branch.lIE	), [1 3 2]);
			VarNxN.Gen.Dfig.Branch.lIF	 = permute(full(	VarM.Gen.Dfig.Branch.lIF	), [1 3 2]);
			VarNxN.Gen.Dfig.Branch.lOR	 = permute(full(	VarM.Gen.Dfig.Branch.lOR	), [1 3 2]);

			VarNxN.Gen.Dfig.Bus.pCF	 = permute(full(	VarM.Gen.Dfig.Bus.pCF	), [1 3 2]);
			VarNxN.Gen.Dfig.Bus.qCF	 = permute(full(	VarM.Gen.Dfig.Bus.qCF	), [1 3 2]);

			VarNxN.Gen.Dfig.Bus.pgE	 = permute(full(	VarM.Gen.Dfig.Bus.pgE	), [1 3 2]);
			VarNxN.Gen.Dfig.Bus.qgE	 = permute(full(	VarM.Gen.Dfig.Bus.qgE	), [1 3 2]);

			VarNxN.Gen.Dfig.Bus.pgR	 = permute(full(	VarM.Gen.Dfig.Bus.pgR	), [1 3 2]);
			VarNxN.Gen.Dfig.Bus.qgR	 = permute(full(	VarM.Gen.Dfig.Bus.qgR	), [1 3 2]);

			VarNxN.Gen.Dfig.Bus.sF	 = permute(full(	VarM.Gen.Dfig.Bus.sF	), [1 3 2]);
			VarNxN.Gen.Dfig.Bus.sR	 = permute(full(	VarM.Gen.Dfig.Bus.sR	), [1 3 2]);

			VarNxN.Gen.Dfig.Bus.xiF	 = permute(full(	VarM.Gen.Dfig.Bus.xiF	), [1 3 2]);
			VarNxN.Gen.Dfig.Bus.xiR	 = permute(full(	VarM.Gen.Dfig.Bus.xiR	), [1 3 2]);

			VarNxN.Gen.Dfig.Bus.vI	 = permute(full(	VarM.Gen.Dfig.Bus.vI	), [1 3 2]);
			VarNxN.Gen.Dfig.Bus.vE	 = permute(full(	VarM.Gen.Dfig.Bus.vE	), [1 3 2]);
			VarNxN.Gen.Dfig.Bus.vF	 = permute(full(	VarM.Gen.Dfig.Bus.vF	), [1 3 2]);
			VarNxN.Gen.Dfig.Bus.vO	 = permute(full(	VarM.Gen.Dfig.Bus.vO	), [1 3 2]);
			VarNxN.Gen.Dfig.Bus.vR	 = permute(full(	VarM.Gen.Dfig.Bus.vR	), [1 3 2]);

			VarNxN.Gen.Dfig.Bus.n_Wnd	 = permute(	VarM.Gen.Dfig.Bus.n_Wnd	, [1 3 2]);
			VarNxN.Gen.Dfig.Bus.P_mecWnd	 = permute(	VarM.Gen.Dfig.Bus.P_mecWnd	, [1 3 2]);
		end
	end

	% Fotovoltaico
	if isfield(VarM, 'Gen')
		if isfield(VarM.Gen, 'Pv')
			VarNxN.Gen.Pv.pPv	 = permute(full(	VarM.Gen.Pv.pPv	), [1 3 2]);
			VarNxN.Gen.Pv.qPv	 = permute(full(	VarM.Gen.Pv.qPv	), [1 3 2]);
			VarNxN.Gen.Pv.s	 = permute(full(	VarM.Gen.Pv.s	), [1 3 2]);
			VarNxN.Gen.Pv.xi	 = permute(full(	VarM.Gen.Pv.xi	), [1 3 2]);
		end
	end

end