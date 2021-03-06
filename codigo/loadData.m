function [Data] = loadData(DflowD)

    Data = [];

    if ~strcmp(DflowD.inFilename,'') && ~strcmp(DflowD.fileCurvaCarga,'')
        [n,et,app] = DimensionDefault(DflowD.inFilename,'bus_data_CVX',DflowD.fileCurvaCarga,DflowD.App);
        Data = DataDefault(n,et,app);

        %% Carga de datos
        % Red
        [Data] = load_distflow_case(Data, DflowD.inFilename, 'bus_data_CVX', 'branch_data_CVX', ...
            DflowD.Trafos, DflowD.Caps, DflowD.Cargas, DflowD.App, DflowD.Switches);

        [Data] = loadCargaCuartHoraria(Data, DflowD.fileCurvaCarga);

        % Tras
        [Data] = cargaTrasDefault(Data, DflowD.Tras);

        % Eolicos
        [Data] = cargaEolicosDefault(Data, DflowD.Eolicos);

        % Fotovoltaicos
        [Data] = cargaPvDefault(Data, DflowD.Solares);

        % Aire Acondicionado
        [Data] = cargaACDefault(Data, DflowD.fileTemp, DflowD.ACs);

        % Baterias
        [Data] = cargaBatDefault(Data, DflowD.Baterias);
        
        % Generador Basico
        [Data] = cargaGBasicoDefault(Data, DflowD.GenBas);
        
        % Utilidad
        [Data] = cargaUtilDefault(Data, DflowD.fileUtilBetaT, DflowD.utilOptFuncCuad, DflowD.Cargas, DflowD.App);

        % Costos
        [Data] = cargaCostosDefault(Data, DflowD.Trafos, DflowD.Caps, DflowD.Switches, DflowD.GenBas, DflowD.fileCostosTension, DflowD.Tras, DflowD.Solares, DflowD.Eolicos);
    end

end