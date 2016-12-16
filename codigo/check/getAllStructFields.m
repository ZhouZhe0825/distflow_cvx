function [gasf, asf] = getAllStructFields(A)

iniStructName = inputname(1);
iniVarStructName = 'A';
currVarStructName = iniVarStructName;
gasf = {currVarStructName};

i = 1;

while i <= length(gasf)
% while ~isempty(gasf)
    currVarStructName = gasf{i};
    eval(['Aux = ', currVarStructName, ';']);
    if isstruct(Aux)
        ind = setdiff(1:length(gasf), [i]);
        gasf = gasf(ind);
        fn = fieldnames(Aux);
        fnComplete = cellfun(@(x)strcat([currVarStructName,'.'],x), fn, 'UniformOutput', false);
        gasf = [gasf;fnComplete];
    else
        i = i+1;
    end

%     currVarStructName1 = gasf{1};
%     eval(['Aux = ', currVarStructName1, ';']);
%     gasf = gasf(2:end);
%     fn1 = fieldnames(Aux);
%     fnComplete = cellfun(@(x)strcat([currVarStructName1,'.'],x), fn1, 'UniformOutput', false);
%     gasf = [gasf;fnComplete];
end
asf = sort(cellfun(@(x)strrep(x, [iniVarStructName,'.'], ''), gasf, 'UniformOutput', false));
gasf = sort(cellfun(@(x)strrep(x, [iniVarStructName,'.'], [iniStructName, '.']), gasf, 'UniformOutput', false));
