function [OutBr] = OutBrMat(T)
	[rowT, colT, ~] = find(T == 1);
	indT = sub2ind(size(T), rowT, colT);
	m = length(indT);
    n = size(T,1);

	OutBr = zeros(n,m);
	for i =1:n
		k = find(T(i,:) == 1)';
		indTk = sub2ind(size(T), repmat(i, [length(k) 1]), k);
		[~,indk,~] = intersect(indT, indTk);
		OutBrI = zeros(1,m);
		OutBrI(indk) = 1;
 		OutBr(i,:) = OutBrI;	
	end
end