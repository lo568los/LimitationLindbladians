% Code to minimize tau for given parameters. This serves as a numerical
% check of the possiblity of thermalization..

beta_list = [1.0,1.2,1.5,1.7,2.0,2.2,2.4,2.7,3.0,3.5];
optimal_value = [];
%load(sprintf("./data_plotting_vsbeta_4/ness_data_NL1=1,NL2=1,NM=2,e=0,beta_r=1,g=0.0016_1.mat"))
%for index = 1:length(beta_list2)
for index1 = 1:12
    load(sprintf("ness_data_NL1 = 1_NL2 = 1_NM = 2_%d.mat",index1))

%decide paramteres.

    NL1 = 1; % has to be atleast 1
    NL2 = 1;
    NM = 2; % has to be atleast 1
    N = NL1+NL2+NM;
    
    dL1 = 2^NL1;
    dL2 = 2^NL2;
    dM = 2^NM;
    d = dL1*dL2*dM;
    e = 0.01;
    
    g = 0.016;
    w0list = zeros(N,1) + 2;
    glist = zeros(N-1,1)- g;
    w0list(4) = 2 + 2*e;
    w0list(5) = 2 + 2*e;
    w0list(6) = 2 + 2*e;
    
    deltalist = zeros(N-1,1);
    beta1 = 0.5;
    ham_type=2;
    %beta2 = beta2;
    
    
    % create_hamiltonian
    H_S = create_hamiltonian(w0list,glist,deltalist,N);
    
    
    
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

    solution = minimize_tauopt_function(NL1,NL2,NM,H_S,dm_ness);
    optimal_value(index) = solution.optimal_val;
    disp(solution.cvx_status);
    
    %rho_th = dm_ness;  %setting our non-eq setup rho_th
    
    %cvx_begin sdp 
     %   cvx_precision high %set CVX precision
     %   variable H_LS1(dL1,dL1) hermitian %declare H_LS to be a hermitian matrix.
     %   variable H_LS2(dL2,dL2) hermitian %declare H_LS to be a hermitian matrix.
     %   variable gamma_matrix1(dL1^2-1,dL1^2-1) hermitian semidefinite %declare Gamma1 to be PSD matrix.
     %   variable gamma_matrix2(dL2^2-1,dL2^2-1) hermitian semidefinite %declare Gamma2 to be PSD matrix.
    
     %   objfunc = 0;
    
     %   for i = 1:d
     %          objfunc = objfunc+ abs(V(:,i)'*(create_L2(rho_th,H_S,H_LS1,gamma_matrix1,F1,NL1,NM + NL2) + create_L2_2(rho_th,H_S,H_LS2,gamma_matrix2,F2,NL2,NM + NL1))*V(:,i)) ;
     %   end
     %   minimize objfunc %minimize objfunc 
    
     %   subject to %constraints..
     %       trace(gamma_matrix1) == 1;
     %       trace(gamma_matrix2) == 1;
            
    %cvx_end
    
    %we compute and save the individual values that are summed up to obtain
    %tau. 
    %diag_values = zeros(d,1);
    %for i = 1:d
    %    diag_values(i) = V(:,i)'*(create_L2(rho_th,H_S,H_LS1,gamma_matrix1,F1,NL1,NM + NL2) + create_L2(rho_th,H_S,H_LS2,gamma_matrix2,F2,NL2,NM + NL1))*V(:,i);
    %end
    
    %save(sprintf("./data_plotting_vsbeta_3/thermal_data2_%d",index1))
    
    %gamma_matrix_approx = round(gamma_matrix,5); % round to 5 places after decimal.
    %xlswrite("gamma_matrix.xlx",gamma_matrix_approx); % can be used to
    %conveniently print out the matrix ..
end
save(sprintf("./data_plotting_vsbeta_4/diag_data_NL1=%d,e=%.2f,beta_r=%.1f,g=%.3f,ham_type=%d.mat",NL1,e,beta1,g,ham_type))