function [utilCarg] = utilidadCarga()

tarifa = [ ...
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    4.5; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
    1.75; ... 
];

social = [ ...
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	0.5; ... 
	4; ... 
	4; ... 
	4; ... 
	4; ... 
];

% utilCarg = .5*tarifa + .5*social;
utilCarg = .5*social;