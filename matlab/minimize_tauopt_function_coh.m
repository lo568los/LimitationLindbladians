% Code to minimize tau for given parameters. This serves as a numerical
% check of the possiblity of thermalization..

beta_list = [];
optimal_value = [];
for index1 = 1:10
    load(sprintf("coh_data_NL = 2,NM = 3_g%d.mat",index1))
    
    %decide paramteres.
    
    NL = 2; % has to be atleast 1
    NM = 3; % has to be atleast 1
    %NL2 = 2;
    N = NL + NM;
    
    dL = 2^NL;
    %dL2 = 2^NL2;
    dM = 2^NM;
    d = 2^N;
    
    
    w0list = zeros(N,1) + 1;
    glist = zeros(N-1,1)+ g;
    
    deltalist = zeros(N-1,1) +1;
    beta_list(index) = beta;
    %beta2 = 1;
    
    
    % create_hamiltonian
    H_S = create_hamiltonian(w0list,glist,deltalist,N);
    
    
    
    [V_unsorted,D_unsorted] = eig(H_S); %not ordered. We order them so that they match the ordering
    %obtained from Python. This is not required, but can be useful for checking
    %both MATLAB and Python Code).
    [temp,ind] = sort(diag(D_unsorted));
    V = V_unsorted(:,ind);
    
    F1 = generate_orthonormal_basis(NL);
    %F2 = generate_orthonormal_basis(NL2);
    
    if (basis_is_orthonormal(F1) == false)
        warning('Basis NOT orthonormal. Something wrong \n');
    end
    
    %if (basis_is_orthonormal(F2) == false)
    %    warning('Basis NOT orthonormal. Something wrong \n');
    %end
    
    
    length_F1 = length(F1); % should be DL^2-1
    %length_F2 = length(F2);
    
    for index = 1:length_F1
        F1{index} = kron(F1{index},eye(dM)/sqrt(dM));
    end
    
    %for index = 1:length_F2
    %    F2{index} = kron(eye(dM*dL1)/sqrt(dM*dL1), F2{index});
    %end
    
    
    if (basis_is_orthonormal(F1) == false)
        warning('Basis NOT orthonormal. Something wrong here! \n');
    end
    
    %if (basis_is_orthonormal(F2) == false)
    %    warning('Basis NOT orthonormal. Something wrong here! \n');
    %end
    
    
    %% Starting the SDP.
    
    rho_th = dm_ness;  %setting our non-eq setup rho_th
    
    cvx_begin sdp 
        cvx_precision high %set CVX precision
        variable H_LS1(dL,dL) hermitian %declare H_LS to be a hermitian matrix.
        %variable H_LS2(dL2,dL2) hermitian %declare H_LS to be a hermitian matrix.
        variable gamma_matrix1(dL^2-1,dL^2-1) hermitian semidefinite %declare Gamma1 to be PSD matrix.
        %variable gamma_matrix2(dL2^2-1,dL2^2-1) hermitian semidefinite %declare Gamma2 to be PSD matrix.
    
    
        objfunc = norm(create_L2(rho_th,H_S,H_LS1,gamma_matrix1,F1,NL,NM) - L2_red);  %norm
    
        %for i = 1:d
         %      objfunc = objfunc+ abs(V(:,i)'*create_L2(rho_th,H_S,H_LS,gamma_matrix,F,NL,NM)*V(:,i)) ;
        %end
        minimize objfunc %minimize objfunc 
    
        subject to %constraints..
            trace(gamma_matrix1) == 1;
            %trace(gamma_matrix2) == 1;
            
    cvx_end

    save(sprintf("./data_plotting_vsg_3/coh_data4_%d",index1))
    
    %we compute and save the individual values that are summed up to obtain
    %tau. 
    %diag_values = zeros(d,1);
    %for i = 1:d
    %    diag_values(i) = V(:,i)'*create_L2(rho_th,H_S,H_LS,gamma_matrix,F,NL,NM)*V(:,i);
    %end
    
    
end


%gamma_matrix_approx = round(gamma_matrix,5); % round to 5 places after decimal.
%xlswrite("gamma_matrix.xlx",gamma_matrix_approx); % can be used to
%conveniently print out the matrix ..