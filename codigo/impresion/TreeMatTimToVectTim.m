function [VT] = TreeMatTimToVectTim(AT, Et, Tree)

	[rowT, colT, ~] = find(Tree == 1);
	n = size(Tree,1);
    
	m = length(rowT);

	indM = (rowT - 1)*m + (1:m)';
	indW = ((1:m)' -1)*n + colT;
	
	M = zeros(m,n);
	W = M';
	M(indM) = 1;
	W(indW) = 1;
	m = size(M,1);
    VT = zeros(m,Et);
    if Et > 1
        for i = 1:Et
            VT(:,i) = (M*AT(:,:,i)*W).*eye(m)*ones(m,1);
        end
    else
        VT(:,1) = (M*AT(:,:)*W).*eye(m)*ones(m,1);
    end
	
end
