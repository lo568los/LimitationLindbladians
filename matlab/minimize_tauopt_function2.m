% Funcion that minimized tau for given input parameters.

function [output] = minimize_tauopt_function2(NL1,NL2,NM,H_S,dm_ness1,L2_red1)
    
   
    N = NL1+NM + NL2;
    dL1 = 2^NL1;
    dL2 = 2^NL2;
    dM = 2^NM;
    d = dL1*dL2*dM;
    
    
    
    [V_unsorted,D_unsorted] = eig(H_S); %not ordered. We order them so that they match the ordering
    %obtained from Python. This is not required, but can be useful for checking
    %both MATLAB and Python Code).
    [temp,ind] = sort(diag(D_unsorted));
    V = V_unsorted(:,ind);
    
    F1 = generate_orthonormal_basis(NL1);
    F2 = generate_orthonormal_basis(NL2);
    
    if (basis_is_orthonormal(F1) == false)
        warning('Basis NOT orthonormal. Something wrong \n');
    end
    
    if (basis_is_orthonormal(F2) == false)
        warning('Basis NOT orthonormal. Something wrong \n');
    end
    
    
    length_F1 = length(F1); % should be DL^2-1
    
    length_F2 = length(F2);
    
    for index = 1:length_F1
        F1{index} = kron(F1{index},eye(dM*dL2)/sqrt(dM*dL2));
    end
    
    for index = 1:length_F2
        F2{index} = kron(eye(dM*dL1)/sqrt(dM*dL1), F2{index});
    end
    
    
    if (basis_is_orthonormal(F1) == false)
        warning('Basis NOT orthonormal. Something wrong here! \n');
    end
    
    if (basis_is_orthonormal(F2) == false)
        warning('Basis NOT orthonormal. Something wrong here! \n');
    end
    
    
    %% Starting the SDP.
    
    rho_th = dm_ness1;  %setting our non-eq setup rho_th
    output = struct(); % Pre-initialize output
    
    
    cvx_begin sdp 
        cvx_precision high %set CVX precision
        variable H_LS1(dL1,dL1) hermitian %declare H_LS to be a hermitian matrix.
        variable H_LS2(dL2,dL2) hermitian %declare H_LS to be a hermitian matrix.
        variable gamma_matrix1(dL1^2-1,dL1^2-1) hermitian semidefinite %declare Gamma1 to be PSD matrix.
        variable gamma_matrix2(dL2^2-1,dL2^2-1) hermitian semidefinite %declare Gamma2 to be PSD matrix.
    
    
        objfunc = norm(create_L2(rho_th,H_S,H_LS1,gamma_matrix1,F1,NL1,NM + NL2) + create_L2_2(rho_th,H_S,H_LS2,gamma_matrix2,F2,NL2,NM + NL1) - L2_red1);  %norm + create_L2_2(rho_th,H_S,H_LS2,gamma_matrix2,F2,NL2,NM + NL1)
    
        %for i = 1:d
         %      objfunc = objfunc+ abs(V(:,i)'*create_L2(rho_th,H_S,H_LS,gamma_matrix,F,NL,NM)*V(:,i)) ;
        %end
        minimize objfunc %minimize objfunc 
    
        subject to %constraints..
            trace(gamma_matrix1) == 1;
            trace(gamma_matrix2) == 1;
            
    cvx_end
        
    
    %Store all relevant info in output..
    output.optimal_val = cvx_optval; %stores tau_opt
    output.gamma_matrix1 = gamma_matrix1; 
    output.gamma_matrix2 = gamma_matrix2;
    %output.diag_values = diag_values;
    output.H_LS1 = H_LS1;
    output.H_LS2 = H_LS2;
    output.cvx_status = cvx_status;
    %output.F = F;
   
    return
end


