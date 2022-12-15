%% Script that computes and plots tauopt as a function of g.

betalist = logspace(log10(0.1),log10(10),10);

NL = 2;
NM = 5;
N = NL+NM;

w0list = zeros(N,1)+1;
deltalist = zeros(N-1,1)+1;
glist = zeros(N-1,1)+0.1;
beta = 1;

debugInfo = {};
solver_status = {};
optimal_value = [];

for index = 1:length(betalist)
    disp(cat(2,'Iteration number :',num2str(index)));
    beta = betalist(index);
    
    solution  =  minimize_tauopt_function(NL,NM,w0list,glist,deltalist,beta);
    debugInfo{index} = solution;
    solver_status{index} = solution.cvx_status;
    optimal_value(index) = solution.optimal_val;
    disp(solution.cvx_status);
end

plot(betalist, optimal_value,'LineWidth',1,'Marker','.','Color','#0072BD','LineStyle','-');
xlabel("$\beta$",'Interpreter','latex');
ylabel("$\tau_{opt}$",'Interpreter','latex')
fontsize(gca,36,"pixels")



% save both data and fig. 
save(cat(2,'./data_plotting_vsbeta_1/data_',num2str(NL),'_',num2str(NM),'.mat'));
savefig( cat(2,'./data_plotting_vsbeta_1/data_',num2str(NL),'_',num2str(NM),'.fig'));

%saveas(gcf,cat(2,'plot',num2str(data_index),'.png'));



