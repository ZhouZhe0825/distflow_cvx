function [AoT] = matOverTime(A)
	Amod = A;
	if size(Amod,3) > 1
		Amod = sum(abs(Amod),3);
	else
		Amod = abs(Amod);
	end
	AoT = sign(Amod);
end