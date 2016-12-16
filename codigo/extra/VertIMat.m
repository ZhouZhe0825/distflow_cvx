function [VertI] = VertIMat(T)
	[rowT, ~, ~] = find(T == 1);
	m = length(rowT);

	VertI = zeros(m,size(T,1));
	for i =1:m
        VertI(i,rowT(i)) = 1;
	end
end