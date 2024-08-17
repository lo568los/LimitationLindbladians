%% Script that computes and plots tauopt as a function of g. This is for hotter bath

th_list = linspace(0.001,8,15);
betalist = [];

for index = 1:length(th_list)
    betalist(index) = 1/th_list(index);
end

NL = 2;
NM = 2;
N = NL+NM;

w0list = zeros(N,1)+1;
%deltalist = zeros(N-1,1)+1;
glist = zeros(N-1,1)+0.0016;


debugInfo = {};
solver_status = {};
optimal_value = [];
gamma_matrices = {};
lamb_shift_matrices = {};

for index = 1:length(betalist)
    disp(cat(2,'Iteration number :',num2str(index))); %Prints out iteration number
    beta = betalist(index);
    
    solution  =  minimize_tauopt_function2(NL,NM,w0list,glist,beta);  %finds the actual solution
    debugInfo{index} = solution;
    solver_status{index} = solution.cvx_status;
    optimal_value(index) = solution.optimal_val;
    gamma_matrices{index} = solution.gamma_matrix;
    lamb_shift_matrices{index} = solution.H_LS;
    disp(solution.cvx_status); %displays whether if solved or some error
end

plot(betalist, optimal_value,'LineWidth',1,'Marker','.','Color','#0072BD','LineStyle','-');
xlabel("$\beta$",'Interpreter','latex');
ylabel("$\tau_{opt}$",'Interpreter','latex')
fontsize(gca,36,"pixels")



% save both data and fig. 
save(cat(2,'./data_plotting_vsbeta_1/data_',num2str(NL),'_',num2str(NM),'_6.mat'));
savefig( cat(2,'./data_plotting_vsbeta_1/data_',num2str(NL),'_',num2str(NM),'_6.fig'));

%saveas(gcf,cat(2,'plot',num2str(data_index),'.png'));



