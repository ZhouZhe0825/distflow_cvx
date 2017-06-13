function [S] = loadCsvDef(filename, fields)
	A = importdata(filename);
	S = cell2struct(num2cell(A.data),A.textdata,1);
    i = 1;
    allChecked = true;
    while i <= length(fields) && allChecked
        allChecked = isfield(S,fields{i});
        i = i + 1;
    end
    if ~allChecked
        S = [];
    end
end