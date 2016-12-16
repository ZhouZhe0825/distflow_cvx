function [Anxnxt] = MxT2NxNxT(VertI, VertJ, Amxt)

	n = size(VertI,2);
    m = size(VertI,1);
    T = size(Amxt,2);
    Anxnxt = zeros(n,n,T);
    for t=1:T
        Anxnxt(:,:,t) = VertI' * (ones(m,1) * Amxt(:,t)' .*eye(m)) *VertJ;
    end

end