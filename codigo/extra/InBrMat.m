function [InBr] = InBrMat(T)
	[rowT, colT, ~] = find(T == 1);
	indT = sub2ind(size(T), rowT, colT);
	m = length(indT);
    n = size(T,1);

	InBr = zeros(n,m);
	for i =1:n
		k = find(T(:,i) == 1)';
		indTk = sub2ind(size(T), k, repmat(i, [1 length(k)]));
		[~,indk,~] = intersect(indT, indTk);
		InBrI = zeros(1,m);
		InBrI(indk) = 1;
 		InBr(i,:) = InBrI;	
	end
end