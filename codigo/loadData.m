function [Data] = loadData(Trafos, Caps, Switches, App, Eolicos, Solares, Baterias, Cargas, ACs, inFilename, fileCurvaCarga, fileUtilBetaT, fileTemp, fileCostosTension, fileCostosTras)

	%% Carga de datos
	% Red
	[Data] = load_distflow_case(inFilename, 'bus_data_CVX', 'branch_data_CVX', Trafos, Caps, Cargas, App, Switches);

	[Data] = loadCargaCuartHoraria(Data, fileCurvaCarga);

	% Eolicos

	[Data] = cargaEolicosDefault(Data, Eolicos);

	% Fotovoltaicos

	[Data] = cargaPvDefault(Data, Solares);

	% Utilidad

	[Data] = cargaUtilDefault(Data, fileUtilBetaT, Cargas, App);

	% Aire Acondicionado

	[Data] = cargaACDefault(Data, fileTemp, ACs);

	% Parametros de Baterias

	% Baterias
	[Data] = cargaBatDefault(Data, Baterias);

	% Costos
	[Data] = cargaCostosDefault(Data, Trafos, Caps, Switches, fileCostosTension, fileCostosTras, Solares, Eolicos);

end