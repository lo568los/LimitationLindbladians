%% Script that computes and plots tauopt as a function of g.

fixed_glist = [0.01;0.025;0.05;0.075;0.100;0.200;0.300]; %The fixed value of g at all sites. 

NL = 2;
NM = 4;
N = NL+NM;

w0list = zeros(N,1)+1;
deltalist = zeros(N-1,1)+1;

beta = 1;

debugInfo = {};
solver_status = {};
optimal_value = [];

for index = 1:length(fixed_glist)
    disp(cat(2,'Iteration number :',num2str(index)));
    g = fixed_glist(index);
    glist = zeros(N-1,1)+g;

    solution  =  minimize_thermal_error_function(NL,NM,w0list,glist,deltalist,beta);
    debugInfo{index} = solution;
    solver_status{index} = solution.cvx_status;
    optimal_value(index) = solution.optimal_val;
    disp(solution.cvx_status);
end

plot(fixed_glist, optimal_value,'LineWidth',1,'Marker','.','Color','#0072BD','LineStyle','-');
xlabel("g");
ylabel("\tau_{opt}")
fontsize(gca,36,"pixels")



% save both data and fig. 
save(cat(2,'./data_plotting_2/data_',num2str(NL),'_',num2str(NM),'.mat'));
savefig( cat(2,'./data_plotting_2/data_',num2str(NL),'_',num2str(NM),'.fig'));

%saveas(gcf,cat(2,'plot',num2str(data_index),'.png'));



