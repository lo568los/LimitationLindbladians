% Code to minimize tau for given parameters. This serves as a numerical
% check of the possiblity of thermalization..

%decide paramteres.

NL = 2; % has to be atleast 1
NM = 2; % has to be atleast 1
N = NL+NM;

dL = 2^NL;
dM = 2^NM;
d = dL*dM;


w0list = zeros(N,1) + 1;
glist = zeros(N-1,1)+0.0016;
%deltalist = zeros(N-1,1) +1;
beta =1;


% create_hamiltonian
H_S = create_hamiltonian2(w0list,glist,N);



[V_unsorted,D_unsorted] = eig(H_S); %not ordered. We order them so that they match the ordering
%obtained from Python. This is not required, but can be useful for checking
%both MATLAB and Python Code).
[temp,ind] = sort(diag(D_unsorted));
V = V_unsorted(:,ind);

F = generate_orthonormal_basis(NL);

if (basis_is_orthonormal(F) == false)
    warning('Basis NOT orthonormal. Something wrong \n');
end


length_F = length(F); % should be DL^2-1

for index = 1:length_F
    F{index} = kron(F{index}, eye(dM)/sqrt(dM));
end


if (basis_is_orthonormal(F) == false)
    warning('Basis NOT orthonormal. Something wrong here! \n');
end


%% Starting the SDP.

rho_th = expm(-beta*H_S) / trace(expm(-beta*H_S));

cvx_begin sdp 
    cvx_precision high %set CVX precision
    variable H_LS(dL,dL) hermitian %declare H_LS to be a hermitian matrix.
    variable gamma_matrix(dL^2-1,dL^2-1) hermitian semidefinite %declare Gamma to be PSD matrix.

    objfunc = 0;

    for i = 1:d
           objfunc = objfunc+ abs(V(:,i)'*create_L2(rho_th,H_S,H_LS,gamma_matrix,F,NL,NM)*V(:,i)) ;
    end
    minimize objfunc %minimize objfunc 

    subject to %constraints..
        trace(gamma_matrix) == 1;
        
cvx_end

%we compute and save the individual values that are summed up to obtain
%tau. 
diag_values = zeros(d,1);
for i = 1:d
    diag_values(i) = V(:,i)'*create_L2(rho_th,H_S,H_LS,gamma_matrix,F,NL,NM)*V(:,i);
end

save("thermal_data1")

%gamma_matrix_approx = round(gamma_matrix,5); % round to 5 places after decimal.
%xlswrite("gamma_matrix.xlx",gamma_matrix_approx); % can be used to
%conveniently print out the matrix ..


