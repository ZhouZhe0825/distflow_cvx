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
	VarNxN.Red.Bus.Tap	 = permute(full(	VarM.Red.Bus.Tap	), [1 3 2]);

	VarNxN.Red.Bus.pC	 = permute(full(	VarM.Red.Bus.pC	), [1 3 2]);
	VarNxN.Red.Bus.qC	 = permute(full(	VarM.Red.Bus.qC	), [1 3 2]);

	VarNxN.Red.Bus.pN	 = permute(full(	VarM.Red.Bus.pN	), [1 3 2]);
	VarNxN.Red.Bus.qN	 = permute(full(	VarM.Red.Bus.qN	), [1 3 2]);

	VarNxN.Red.Bus.pG	 = permute(full(	VarM.Red.Bus.pG	), [1 3 2]);
	VarNxN.Red.Bus.qG	 = permute(full(	VarM.Red.Bus.qG	), [1 3 2]);

	VarNxN.Red.Bus.qCp	 = permute(full(	VarM.Red.Bus.qCp	), [1 3 2]);
	VarNxN.Red.Bus.Cap	 = permute(full(	VarM.Red.Bus.Cap	), [1 3 2]);

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

			VarNxN.Gen.Dfig.Branch.P	 = zeros(5,5, Etapas,lenWN);
			VarNxN.Gen.Dfig.Branch.P(1,2,:,:)	 = 	VarM.Gen.Dfig.Branch.PIE	;
			VarNxN.Gen.Dfig.Branch.P(1,3,:,:)	 = 	VarM.Gen.Dfig.Branch.PIF	;
			VarNxN.Gen.Dfig.Branch.P(4,5,:,:)	 = 	VarM.Gen.Dfig.Branch.POR	;

			VarNxN.Gen.Dfig.Branch.Q	 = zeros(5,5, Etapas,lenWN);
			VarNxN.Gen.Dfig.Branch.Q(1,2,:,:)	 = 	VarM.Gen.Dfig.Branch.QIE	;
			VarNxN.Gen.Dfig.Branch.Q(1,3,:,:)	 = 	VarM.Gen.Dfig.Branch.QIF	;
			VarNxN.Gen.Dfig.Branch.Q(4,5,:,:)	 = 	VarM.Gen.Dfig.Branch.QOR	;

			VarNxN.Gen.Dfig.Branch.l	 = zeros(5,5, Etapas,lenWN);
			VarNxN.Gen.Dfig.Branch.l(1,2,:,:)	 = 	VarM.Gen.Dfig.Branch.lIE	;
			VarNxN.Gen.Dfig.Branch.l(1,3,:,:)	 = 	VarM.Gen.Dfig.Branch.lIF	;
			VarNxN.Gen.Dfig.Branch.l(4,5,:,:)	 = 	VarM.Gen.Dfig.Branch.lOR	;

			VarNxN.Gen.Dfig.Bus.v	 = zeros(5,1, Etapas,lenWN);
			VarNxN.Gen.Dfig.Bus.v(1,1,:,:)	 = 	VarM.Gen.Dfig.Bus.vI	;
			VarNxN.Gen.Dfig.Bus.v(2,1,:,:)	 = 	VarM.Gen.Dfig.Bus.vE	;
			VarNxN.Gen.Dfig.Bus.v(3,1,:,:)	 = 	VarM.Gen.Dfig.Bus.vF	;
			VarNxN.Gen.Dfig.Bus.v(4,1,:,:)	 = 	VarM.Gen.Dfig.Bus.vO	;
			VarNxN.Gen.Dfig.Bus.v(5,1,:,:)	 = 	VarM.Gen.Dfig.Bus.vR	;

			VarNxN.Gen.Dfig.Bus.pC	 = zeros(5,1, Etapas,lenWN);
			VarNxN.Gen.Dfig.Bus.pC(3,1,:,:)	 = 	VarM.Gen.Dfig.Bus.pCF	;

			VarNxN.Gen.Dfig.Bus.qC	 = zeros(5,1, Etapas,lenWN);
			VarNxN.Gen.Dfig.Bus.qC(3,1,:,:)	 = 	VarM.Gen.Dfig.Bus.qCF	;

			VarNxN.Gen.Dfig.Bus.pg	 = zeros(5,1, Etapas,lenWN);
			VarNxN.Gen.Dfig.Bus.pg(2,1,:,:)	 = 	VarM.Gen.Dfig.Bus.pgE	;
			VarNxN.Gen.Dfig.Bus.pg(5,1,:,:)	 = 	VarM.Gen.Dfig.Bus.pgR	;

			VarNxN.Gen.Dfig.Bus.qg	 = zeros(5,1, Etapas,lenWN);
			VarNxN.Gen.Dfig.Bus.qg(2,1,:,:)	 = 	VarM.Gen.Dfig.Bus.qgE	;
			VarNxN.Gen.Dfig.Bus.qg(5,1,:,:)	 = 	VarM.Gen.Dfig.Bus.qgR	;

			VarNxN.Gen.Dfig.Bus.s	 = zeros(5,1, Etapas,lenWN);
			VarNxN.Gen.Dfig.Bus.s(3,1,:,:)	 = 	VarM.Gen.Dfig.Bus.sF	;
			VarNxN.Gen.Dfig.Bus.s(5,1,:,:)	 = 	VarM.Gen.Dfig.Bus.sR	;

			VarNxN.Gen.Dfig.Bus.xi	 = zeros(5,1, Etapas,lenWN)
			VarNxN.Gen.Dfig.Bus.xi(3,1,:,:)	 = 	VarM.Gen.Dfig.Bus.xiF	;
			VarNxN.Gen.Dfig.Bus.xi(5,1,:,:)	 = 	VarM.Gen.Dfig.Bus.xiR	;

			VarNxN.Gen.Dfig.Bus.n_Wnd	 = permute(	VarM.Gen.Dfig.Bus.n_Wnd	, [4 3 2 1]);
			VarNxN.Gen.Dfig.Bus.P_mecWnd	 = permute(	VarM.Gen.Dfig.Bus.P_mecWnd	, [4 3 2 1]);
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