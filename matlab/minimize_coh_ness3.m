% Code to minimize tau for given parameters. This serves as a numerical
% check of the possiblity of thermalization..

%betal_list = [0.100 0.15199111 0.23101297 0.35111917 0.53366992 0.81113083 1.23284674 1.87381742 2.84803587 4.32876128 6.57933225 10.000];
%g_list2 = [0.010,0.019,0.035,0.066,0.123,0.231,0.433,1.52,2.84803587,5.33669923,10.000];
g_list = [0.0010,0.00187382,0.00351119,0.00657933,0.01232847,0.0231013,0.04328761,0.08111308,0.15199111,0.28480359,0.53366992,1.0000];
betar_list = [1.0];
e_list = [0.00];
h_list = [1];
optimal_value = [];
%load(sprintf("./data_plotting_vsbeta_4/ness_data_NL1=1,NL2=1,NM=2,e=0,beta_r=1,g=0.0016_1.mat"))
%for index = 1:length(beta_list2)
%ham_type=1;
%e = 0.00;
for e1 = 1:length(e_list)
    for h1 = 1:length(h_list)
        ham_type = h_list(h1);
        e = e_list(e1);
        for index1 = 1:length(betar_list)
            disp(betar_list(index1))
            for index2 = 1:length(g_list)
                load(sprintf("./Other_coh_data/coh_data_NL = 3,NM = 1_g%d.mat",index2))

                %ness_data_NL1=2,NL2=2,NM=2,e=%.2f,beta_r=%.1f,beta_l=1.0,g=%.4f_%d.mat
                
            
            %decide paramteres.
            
                NL1 = 3; % has to be atleast 1
                NL2 = 0;
                NM = 1; % has to be atleast 1
                N = NL1+NL2+NM;
                
                dL1 = 2^NL1;
                dL2 = 2^NL2;
                dM = 2^NM;
                d = dL1*dL2*dM;
                
                
                g = g_list(index2);
                w0list = zeros(N,1) + 1;
                glist = zeros(N-1,1)+ g;
                w0list(4) = 1 + e;
                w0list(5) = 1 + e;
                w0list(6) = 1 + e;
                if ham_type == 1
                    deltalist = zeros(N-1,1) + 1;
                end
                if ham_type == 2
                    deltalist = zeros(N-1,1);
                end
                %beta1 = betal_list(index1);
                
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
            
                solution = minimize_tauopt_function2(NL1,NL2,NM,H_S,dm_ness,L2_red);
                optimal_value(index2) = solution.optimal_val;
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
            save(sprintf("./data_plotting_vsg_3/coh_data_NL1=%d,NM=0,e=%.2f,beta_r=%.1f,beta_l=1.0,ham_type=%d.mat",NL1,e,betar_list(index1),ham_type))
        end
    end
end