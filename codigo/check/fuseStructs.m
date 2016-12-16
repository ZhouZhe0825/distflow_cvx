function [res] = fuseStructs(Astruct, Bstruct)

[Agasf, ~] = getAllStructFields(Astruct);
[Bgasf, ~] = getAllStructFields(Bstruct);
res = Astruct;
indA = (1:size(Agasf,1));
indB = (1:size(Bgasf,1));

i = 1;
while ~isempty(indB)
	indI = indB(i);
	Bfield = Bgasf{indI};
	j = 1;
	find = false;
	while ~isempty(indA) && j <= length(indA) && ~find
		indJ = indA(j);
		Afield = Agasf{indJ};
		if strcmp(strrep(Bfield, 'Bstruct.', ''), strrep(Afield, 'Astruct.', ''))
			BAux = [];
			AAux = [];
			eval(['BAux = ', Bfield,';']);
			eval(['AAux = ', Afield,';']);
			AAux = AAux+BAux; % TODO: chequear que son de iguales dimensiones 
			eval(['res.', strrep(Bfield, 'Bstruct.', ''), ' = AAux;']);
            find = true;
			indB(i) = [];
			indA(j) = [];
		else
			j = j+1;
		end
	end
	if ~find
		BAux = [];
		eval(['BAux = ', Bfield,';']);
		eval(['res.', strrep(Bfield, 'Bstruct.', ''), ' = BAux;']);
		indB(i) = [];
	end
end
