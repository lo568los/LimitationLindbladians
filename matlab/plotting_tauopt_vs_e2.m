%% Script that computes and plots tauopt as a function of g. This is for colder bath
tc_list = [0.1];
betalist = [];

for index = 1:length(tc_list)
    betalist(index) = 1/tc_list(index);
end


NL = 2;
NM = 4;
N = NL+NM;

w0list = zeros(N,1)+1;
w0list2 = zeros(N,1)+1;
%deltalist = zeros(N-1,1)+1;
glist = zeros(N-1,1)+0.0016;

deltalist = [0, 0.00263158, 0.00526316, 0.00789474, 0.01052632, 0.01315789, 0.01578947, 0.01842105, 0.02105263, 0.02368421, 0.02631579, 0.02894737, 0.03157895, 0.03421053, 0.03684211, 0.03947368, 0.04210526, 0.04473684, 0.04736842, 0.05];


super_optimal_value = [];
super_gamma_matrices = {};
super_lamb_shift_matrices = {};

for index = 1:length(betalist)
    beta = betalist(index);
    optimal_value = [];
    gamma_matrices = {};
    lamb_shift_matrices = {};
    for index2 = 1:length(deltalist)
        disp(cat(2,'Iteration number(Th) :',num2str(index),' Iteration number(e) :',num2str(index2))); %Prints out iteration number
        
        delta = deltalist(index2);
        %output = modified_list(N,w0list,delta);
        w0list2(4) = w0list(4) + delta;
        w0list2(5) = w0list(5) + delta;
        w0list2(6) = w0list(6) + delta;

        disp(w0list2);

        
        solution  =  minimize_tauopt_function2(NL,NM,w0list2,glist,beta);  %finds the actual solution
        optimal_value(index2) = solution.optimal_val;
        gamma_matrices{index2} = solution.gamma_matrix;
        lamb_shift_matrices{index2} = solution.H_LS;
        disp(solution.cvx_status);
        disp(solution.optimal_val);%displays whether if solved or some error
    end
    super_gamma_matrices{index} = gamma_matrices;
    super_lamb_shift_matrices{index} = lamb_shift_matrices;
end

%plot(betalist, optimal_value,'LineWidth',1,'Marker','.','Color','#0072BD','LineStyle','-');
%xlabel("$\beta$",'Interpreter','latex');
%ylabel("$\tau_{opt}$",'Interpreter','latex')
%fontsize(gca,36,"pixels")



% save both data and fig. 
save(cat(2,'./data_plotting_vs_e/data_',num2str(NL),'_',num2str(NM),'_4.mat'));
%savefig( cat(2,'./data_plotting_vsbeta_1/data_',num2str(NL),'_',num2str(NM),'_5.fig'));

%saveas(gcf,cat(2,'plot',num2str(data_index),'.png'));



