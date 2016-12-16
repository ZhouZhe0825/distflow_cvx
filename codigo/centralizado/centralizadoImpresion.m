function [Var, opt, status, Dmod] = centralizadoImpresion(Trafos, Caps, Cargas, outFilename, Data, Config, utilCarg)
                               
	cantTrafos = length(Trafos);
	cantCaps = length(Caps);
	cantCargs = length(Cargas);

	%% Resolucion

	[Var, opt, status, Dmod] = llamarCentralizado(Data, Config, utilCarg);

	printSalidasDistflow(Var, Dmod, Config, cantTrafos, cantCaps, cantCargs, outFilename, [], [], [], [], []);

end