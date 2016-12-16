function [P_mec, n_] = calculoPotenciaEolica(v, vmm, vm, vM, vMM, Omega, G, P_nMec, Np, R_, rho, ws, c1, c2, c3, c4, c5, c6, c7, lambda_opt)

	lambda = 0; %TODO Poner bien
	beta = 0; %TODO Poner bien

    P_mec = zeros(size(v));
    n_ = zeros(size(v));
    
    for i = 1:length(v)
        if v(i) > vmm && v(i) <= vMM
            if v(i) > vmm && v(i) <= vm
                beta = 0;
                lambda = lambda_opt;
            elseif v(i) > vm
                beta = 0;
                lambda = (R_ * Omega) / (v(i) * G);
            end
            n_(i) = 1 - (1/ws)*(Np * v(i) * G * lambda/ R_);
            P_mec(i) = (1/2) * rho * pi * R_^2 * calculoC_p(lambda, beta, c1, c2, c3, c4, c5, c6, c7)*v(i)^3;
            if v(i) <= vMM && v(i) > vM
                P_mec(i) = P_nMec;
            end
        else
            n_(i) = 0;
            P_mec(i) = 0;
        end
        
    end
	
    P_mec = P_mec/1e6;

end

function [val] = calculoC_p(lambda, beta, c1, c2, c3, c4, c5, c6, c7)

	inv_lambda_i = 1/(lambda+c6*beta) - c7/(beta^3+1);
	val = c1*(c2*inv_lambda_i - c3*beta - c4) * exp(-c5*inv_lambda_i);
	
end