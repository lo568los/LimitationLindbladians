function output=create_L2(rho_th,H_S,H_LS,gamma,F,NL,NM)
    dL = 2^NL;
    dM = 2^NM;


    
    Heff = kron(H_LS,eye(dM));
    term1 = 1.0i * (rho_th*Heff -Heff*rho_th );
    term2 = 0;

    for i = 1:(dL^2-1)
        for j = 1: (dL^2-1)
            term2 = term2 + gamma(i,j)* (F{j}*rho_th*F{i}'-0.5*F{i}'*F{j}*rho_th - 0.5*rho_th*F{i}'*F{j} );
        end
    end

    output = term1+term2;
    
end