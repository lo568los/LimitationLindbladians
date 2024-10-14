%% Script that computes and plots tauopt as a function of g. This is for colder bath
th_list = [0.4,0.7,1.0];
betalist = [];

for index = 1:length(th_list)
    betalist(index) = 1/th_list(index);
end


NL = 2;
NM = 4;
N = NL+NM;

delta = 0.01;

w0list = zeros(N,1)+1;
w0list2 = zeros(N,1)+1;
w0list2(4) = w0list(4) + delta;
w0list2(5) = w0list(5) + delta;
w0list2(6) = w0list(6) + delta;
%deltalist = zeros(N-1,1)+1;


gvals = [0.001, 0.00147368, 0.00194737, 0.00242105, 0.00289474, 0.00336842, 0.00384211, 0.00431579, 0.00478947, 0.00526316, 0.00573684, 0.00621053, 0.00668421, 0.00715789, 0.00763158, 0.00810526, 0.00857895, 0.00905263, 0.00952632, 0.01];


super_optimal_value = [];
super_gamma_matrices = {};
super_lamb_shift_matrices = {};

for index = 1:length(betalist)
    beta = betalist(index);
    optimal_value = [];
    gamma_matrices = {};
    lamb_shift_matrices = {};
    for index2 = 1:length(gvals)
        disp(cat(2,'Iteration number(Th) :',num2str(index),' Iteration number(e) :',num2str(index2))); %Prints out iteration number

        glist = zeros(N-1,1)+gvals(index2);
        
        
        %output = modified_list(N,w0list,delta);
        disp(glist)
        
        solution  =  minimize_tauopt_function(NL,NM,w0list2,glist,beta);  %finds the actual solution
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
save(cat(2,'./data_plotting_vs_g/data_',num2str(NL),'_',num2str(NM),'_1.mat'));
%savefig( cat(2,'./data_plotting_vsbeta_1/data_',num2str(NL),'_',num2str(NM),'_5.fig'));

%saveas(gcf,cat(2,'plot',num2str(data_index),'.png'));