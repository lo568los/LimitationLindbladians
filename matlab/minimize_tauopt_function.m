% Funcion that minimized tau for given input parameters.

function [output] = minimize_tauopt_function(NL,NM,w0list,glist,deltalist,beta)
    
   
    N = NL+NM;
    dL = 2^NL;
    dM = 2^NM;
    d = dL*dM;
    
    
    
    H_S = create_hamiltonian(w0list,glist,deltalist,N);
    
    
    [V,~] = eig(H_S);
    
    F = generate_orthonormal_basis(NL);
    
    if (basis_is_orthonormal(F) == false)
        warning('Basis NOT orthonormal. Something wrong \n');
    end
    
    
    length_F = length(F); % should be DL^2-1
    %F = F{1:length_F-1};
    
    for index = 1:length_F
        F{index} = kron(F{index}, eye(dM)/sqrt(dM));
    end
    
    
    if (basis_is_orthonormal(F) == false)
        warning('Basis NOT orthonormal. Something wrong here! \n');
    end
    
    
    %% Starting the SDP.
    
    rho_th = expm(-beta*H_S) / trace(expm(-beta*H_S));
    
    cvx_begin sdp quiet
        cvx_precision high
        variable H_LS(dL,dL) hermitian
        variable gamma_matrix(dL^2-1,dL^2-1) hermitian semidefinite
    
        objfunc = 0;
    
        for i = 1:d
               objfunc = objfunc+ abs(V(:,i)'*create_L2(rho_th,H_S,H_LS,gamma_matrix,F,NL,NM)*V(:,i)) ;
        end
        minimize objfunc
    
        subject to
            trace(gamma_matrix) == 1;
            
    cvx_end
        
    diag_values = zeros(d,1);
    for i = 1:d
        diag_values(i) = V(:,i)'*create_L2(rho_th,H_S,H_LS,gamma_matrix,F,NL,NM)*V(:,i);
    end
    

    %Store all relevant info in output..
    output.optimal_val = cvx_optval; %stores tau_opt
    output.gamma_matrix = gamma_matrix; 
    output.diag_values = diag_values;
    output.H_LS = H_LS;
    output.cvx_status = cvx_status;

    return
end


