%% reads from the .mat file and computes \tau for a single set of modified parameters.

load('thermal_data2.mat');
glist(5) = 0.2;  %modify parameters..



H_S = create_hamiltonian2(w0list,glist,N); % hamiltonian
[V_unsorted,D_unsorted] = eig(H_S); %not ordered . need ordered to check with Python
[temp,ind] = sort(diag(D_unsorted));
V = V_unsorted(:,ind);

rho_th = expm(-beta*H_S) / trace(expm(-beta*H_S));

tau = 0;

for i = 1:d
       tau = tau + abs(V(:,i)'*create_L2_2(rho_th,H_S,H_LS,gamma_matrix,F,NL,NM)*V(:,i)) ;
end

tau