%% code for plotting tau vs g(bulk site) for various configurations..


fixed_glist = logspace(log10(0.001),log10(1),10);


% This loads the Gamma, H_LS, H_S etc.
tau_list = [];
%site_loc = 4; %The location of g that we want to change value at. 
NL = 3;
NM = 3;

%modify other parameters 
%deltalist(3) = 0.4;
%deltalist(4) = 1.2;
%glist(3) = 0.3;

for index = 1:length(fixed_glist)
    load(sprintf('./Other_coh_data/coh_data_NL = %d,NM = %d_g%d.mat',NL,NM,index)); % load the data.
    disp(cat(2,'Iteration number :',num2str(index)));
    g = fixed_glist(index);
    %glist(site_loc) = fixed_glist(index); %Modify g at that location.

    %create new H_S for modified parameters
    H_S = create_hamiltonian(w0list,glist,deltalist,N); % hamiltonian
    [V_unsorted,D_unsorted] = eig(H_S); %not ordered . need ordered to check with Python
    [temp,ind] = sort(diag(D_unsorted));
    V = V_unsorted(:,ind);

    beta = 1;

    % create new rho_th for modified parameters
    rho_th = expm(-beta*H_S) / trace(expm(-beta*H_S));

    %start computing tau (although we call this objfunc, there is no
    %optimization involved here)

    solution = minimize_tauopt_function2(NL,0,NM,H_S,dm_ness,L2_red);
    tau_list(index) = solution.optimal_val;

    %objfunc = 0;

    %for i = 1:d
    %       objfunc = objfunc+ abs(V(:,i)'*create_L2(rho_th,H_S,H_LS,gamma_matrix,F,NL,NM)*V(:,i)) ;
    %end
    
    
            
end
save("./data_plotting_vsg_3/coh_data_NL1=3,NM=3.mat")

%plot(fixed_glist, tau_list,'LineWidth',1,'Marker','.','Color','#0072BD','LineStyle','-');
%xlabel("g");
%ylabel("\tau")
%fontsize(gca,36,"pixels")

%legend()



%save('./tau_plotting/data_5.mat');


%saveas(gcf,cat(2,'plot',num2str(data_index),'.png'));


