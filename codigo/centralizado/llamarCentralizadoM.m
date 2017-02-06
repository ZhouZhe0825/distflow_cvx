function [Var, opt, DataM] = llamarCentralizadoM(Data, Config)


    DataM = Data;
    [DataM] = reshapeDataM(DataM, Config);

	[Var, opt] = distflowCentralizadoM(DataM, Config);

% 	DataM_f = DataM;
% 
% 	leyenda = 'Segunda vuelta centralizado'
% 
% 	DataM_f.Fixed.Cap = round(Var.Red.Bus.Cap);
% 	DataM_f.Fixed.Tap = round(Var.Red.Bus.Tap); 
% 	DataM_f.Fixed.stCh = round(Var.ClNI.start);
% 	DataM_f.Fixed.onCh = round(Var.ClNI.on);
% 	DataM_f.Fixed.y = round(Var.Red.Branch.y);
% 	DataM_f.Fixed.z = round(Var.Red.Branch.z);
% 
% 	[Var, opt, status] = distflow_whole_sw_d(DataM_f, Config, true, util);
% 
% 	leyenda = ['Inicial Duales - ' status]
% 
% 	leyenda = 'Fin centralizado'
% 	save(Config.workspace_var_file);
end