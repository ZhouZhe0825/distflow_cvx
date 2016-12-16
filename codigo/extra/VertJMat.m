function [VertJ] = VertJMat(T)
	[~, colT, ~] = find(T == 1);
	m = length(colT);

	VertJ = zeros(m,size(T,1));
	for i =1:m
        VertJ(i,colT(i)) = 1;
	end
end