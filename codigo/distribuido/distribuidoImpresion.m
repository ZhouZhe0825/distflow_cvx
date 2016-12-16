function [Var_dist_conE, Var_centr, opt_dist_conE, opt_centr, status, Dmod, Ev] = distribuidoImpresion(Trafos, Caps, Cargas, outFilename, Data, Config, utilCarg, Var_centr, opt_centr, Var_ini, opt_ini)

% 	workspace_var_file = [outFilename, '_var.mat'];
	outFilename_centr = [outFilename, '_centr'];
	outFilename_dist_conE = [outFilename, '_dist_conE'];

	cantTrafos = length(Trafos);
	cantCaps = length(Caps);
	cantCargs = length(Cargas);

	% cantCases = length(Cases);

	%% Resolucion

% 	Output = initSalidas_whole_sw_dist(Data, Config, Caps, Cargas, cantTrafos, cantCaps, cantCargs);
% 	Output_centr = initSalidas_whole_sw_dist(Data, Config, Caps, Cargas, cantTrafos, cantCaps, cantCargs);
% 	Output_dist_conE = initSalidas_whole_sw_dist(Data, Config, Caps, Cargas, cantTrafos, cantCaps, cantCargs);
    
    [Var_dist_conE, Var_centr, Var_F, opt_dist_conE, opt_centr, opt_F, status, Dmod, Ev] = llamarDistribuido(Data, Config, utilCarg, Var_centr, opt_centr, Var_ini, opt_ini);
    
%     [Output] =           actualizarSalidas_whole_sw_dist(Var_F,         Dmod, Config, Output);
%     [Output_centr] =     actualizarSalidas_whole_sw_dist(Var_centr,     Dmod, Config, Output_centr);
%     [Output_dist_conE] = actualizarSalidas_whole_sw_dist(Var_dist_conE, Dmod, Config, Output_dist_conE);

	printSalidasDistflow(Var_F,         Dmod, Config, cantTrafos, cantCaps, cantCargs, outFilename,           Ev.opt, Ev.mu, Ev.lambda, Ev.DifP, Ev.DifQ);
	printSalidasDistflow(Var_centr,     Dmod, Config, cantTrafos, cantCaps, cantCargs, outFilename_centr,     [],    [],   [],       [],     []);
	printSalidasDistflow(Var_dist_conE, Dmod, Config, cantTrafos, cantCaps, cantCargs, outFilename_dist_conE, [],    [],   [],       [],     []);
end