function [res] = checkEqualStructs(Astruct, Bstruct, AstructName, BstructName, tol)

[Agasf, ~] = getAllStructFields(Astruct);
[Bgasf, ~] = getAllStructFields(Bstruct);
res = [];
if (size(Agasf,1) == size(Bgasf,1)) && (size(Agasf,2) == size(Bgasf,2)) 
    indA = (1:size(Agasf,1));
    indB = (1:size(Bgasf,1));
    
    i = 1;
    while ~isempty(indA)
        indI = indA(i);
        Afield = Agasf{indI};
        j = 1;
        find = false;
        while ~isempty(indB) && ~find
            indJ = indB(j);
            Bfield = Bgasf{indJ};
            if strcmp(strrep(Afield, 'Astruct.', ''), strrep(Bfield, 'Bstruct.', ''))
                disc = [];
                AAux = [];
                BAux = [];
                eval(['AAux = ', Afield,';']);
                eval(['BAux = ', Bfield,';']);
                disc = similarMat(AAux, BAux, tol);
                if ~isempty(disc)
                    d = mat2cell(disc, ones(size(disc,1),1), ones(1,size(disc,2)));
                    if ~isempty(AstructName)
                        Afield = strrep(Afield,'Astruct', AstructName);
                    end
                    if ~isempty(BstructName)
                        Bfield = strrep(Bfield,'Bstruct', BstructName);
                    end
                    res = [res;{Afield, Bfield, ''};d];
                end
                find = true;
                indA(i) = [];
                indB(j) = [];
            else
                j = j+1;
            end
        end
        if ~find
            res = [res;{Afield, 'none', ''};{Inf Inf Inf}];
            indA(i) = [];
        end
    end
else
    % distinto numero de campos, no son iguales
end

end