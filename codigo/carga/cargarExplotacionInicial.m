function [yIni] = cargarExplotacionInicial(DataM)

	DataIni = DataM;

	DataIni.Cost.cdv = DataM.Cost.cdv(:,1);
	DataIni.Cost.cY = DataM.Cost.cY(:,1);
	DataIni.Cost.delta = DataM.Cost.delta(:,1);
	DataIni.Cost.m = DataM.Cost.m(:,1);
	DataIni.Red.Branch.lTop = DataM.Red.Branch.lTop(:,1);
	DataIni.Red.Branch.NtrLow = DataM.Red.Branch.NtrLow(:,1);
	DataIni.Red.Branch.NtrTop = DataM.Red.Branch.NtrTop(:,1);
	DataIni.Red.Branch.Tap = DataM.Red.Branch.Tap(:,1);
	DataIni.Red.Branch.r = DataM.Red.Branch.r(:,1);
	DataIni.Red.Branch.x = DataM.Red.Branch.x(:,1);
	DataIni.Red.Branch.yLow = DataM.Red.Branch.yLow(:,1);
	DataIni.Red.Branch.yTop = DataM.Red.Branch.yTop(:,1);
	DataIni.Red.Bus.Cap = DataM.Red.Bus.Cap(:,1);
	DataIni.Red.Bus.NcpLow = DataM.Red.Bus.NcpLow(:,1);
	DataIni.Red.Bus.NcpTop = DataM.Red.Bus.NcpTop(:,1);
	DataIni.Red.Bus.uLow = DataM.Red.Bus.uLow(:,1);
	DataIni.Red.Bus.uTop = DataM.Red.Bus.uTop(:,1);
	DataIni.Red.Bus.pCLow = DataM.Red.Bus.pCLow(:,1);
	DataIni.Red.Bus.qCLow = DataM.Red.Bus.qCLow(:,1);
	DataIni.Gen.Tras.pgLow = DataM.Gen.Tras.pgLow(:,1);
	DataIni.Gen.Tras.pgTop = DataM.Gen.Tras.pgTop(:,1);
	DataIni.Gen.Tras.qgLow = DataM.Gen.Tras.qgLow(:,1);
	DataIni.Gen.Tras.qgTop = DataM.Gen.Tras.qgTop(:,1);
	
	leyenda = '------------------------ Calculo de explotacion ------------------------'
	
	[Var, status] = df_InicialM(DataIni);

	leyenda = status

	yIni = Var.Red.Branch.y;



