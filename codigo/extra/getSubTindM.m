function [SubTM, SubTM_T] = getSubTindM(T,SubT)

	[rSubT, cSubT, ~] = find(SubT == 1);

	Tind = find(T == 1);
	SubTind = sub2ind(size(T), rSubT, cSubT);
	SubTM_Tind = sub2ind(size(T), cSubT, rSubT);

	SubTM = zeros(size(SubTind));
	SubTM_T = zeros(size(SubTind));
	for i = 1:length(SubTind)
		SubTM(i) = find(Tind == SubTind(i));
		SubTM_T(i) = find(Tind == SubTM_Tind(i));
	end
end
