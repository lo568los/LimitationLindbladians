% Code to minimize tau for given parameters. This serves as a numerical
% check of the possiblity of thermalization..

NL1 = 1; % has to be atleast 1
NL2 = 1;
NM = 2; % has to be atleast 1
N = NL1+NL2+NM;
e=0;
beta_r=0.5;
g=0.0016;
ham_type=1;

beta_list = [0.1,0.15199111,0.23101297,0.35111917,0.53366992,0.81113083,1.23284674,1.87381742,2.84803587,4.32876128,6.5793322,10.0];
optimal_value = [];
%ness_data_NL1=2,NL2=2,NM=2,e=0.01,beta_r=10,g=0.016_2
load(sprintf("./data_plotting_vsbeta_4/ness_data_NL1=%d,NL2=%d,NM=%d,e=%d,beta_r=%.1f,g=%.4f_%d.mat",NL1,NL2,NM,e,beta_r,g,ham_type))
for index = 1:length(beta_list)
%decide paramteres.

    
    
    dL1 = 2^NL1;
    dL2 = 2^NL2;
    dM = 2^NM;
    d = dL1*dL2*dM;
    
    
    
    w0list = zeros(N,1) + 2;
    glist = zeros(N-1,1)+ -g;
    w0list(4) = 2 + 2*e;
    w0list(5) = 2 + 2*e;
    w0list(6) = 2 + 2*e;
    
    deltalist = zeros(N-1,1) +1;
    beta1 = beta_r;
    beta2 = beta_list(index);
    
    
    % create_hamiltonian
    if ham_type ==1
        deltalist = zeros(N,1) + 1;
        H_S = create_hamiltonian(w0list,glist,deltalist,N);
    end

    if ham_type ==2
        deltalist = zeros(N,1);
        H_S = create_hamiltonian(w0list,glist,deltalist,N);
    end
    
    
    dm_ness1 = dm_ness(index);
    L2_red1 = L2_red(index);
    solution = minimize_tauopt_function2(NL1,NL2,NM,H_S,dm_ness1,L2_red1);
    optimal_value(index) = solution.optimal_val;
    disp(solution.cvx_status);
    
    
    
    %gamma_matrix_approx = round(gamma_matrix,5); % round to 5 places after decimal.
    %xlswrite("gamma_matrix.xlx",gamma_matrix_approx); % can be used to
    %conveniently print out the matrix ..
end

save(sprintf("./data_plotting_vsbeta_4/coh_data_NL1=%d_NL2=%d_NM=%d,e=%f,beta_r=%d,g=%f_%d.mat",NL1,NL2,NM,e,beta_r,g,ham_type))