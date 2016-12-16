function [Amxt] = NxNxT2MxT(VertI, VertJ, Anxnxt)

	m = size(VertI,1);
    T = size(Anxnxt,3);
    Amxt = zeros(m,T);
    for t=1:T
        Amxt(:,t) = (VertI * Anxnxt(:,:,t) * VertJ').*eye(m)*ones(m,1);
    end

end