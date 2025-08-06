e = 0.00;
ham_type=1;
NL1=2;
g=0.0100;


X1=getfield(load(sprintf('./coh_data_new2_NL1=%d,e=%.2f,beta_r=10.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'s_list');
Y1=getfield(load(sprintf('./coh_data_new2_NL1=%d,e=%.2f,beta_r=10.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'optimal_value');

X2=getfield(load(sprintf('./coh_data_new2_NL1=%d,e=%.2f,beta_r=5.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'s_list');
Y2=getfield(load(sprintf('./coh_data_new2_NL1=%d,e=%.2f,beta_r=5.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'optimal_value');

X3=getfield(load(sprintf('./coh_data_new2_NL1=%d,e=%.2f,beta_r=1.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'s_list');
Y3=getfield(load(sprintf('./coh_data_new2_NL1=%d,e=%.2f,beta_r=1.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'optimal_value');

X4=getfield(load(sprintf('./coh_data_new2_NL1=%d,e=%.2f,beta_r=0.5,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'s_list'); 
Y4=getfield(load(sprintf('./coh_data_new2_NL1=%d,e=%.2f,beta_r=0.5,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'optimal_value');

%% make sure all these elists are the same !!! Makes comparison easier!!




semilogy(X1,Y1,'DisplayName',"$\beta_R = 10.0$",'LineWidth',5,'Marker','*','MarkerSize',20,'LineStyle','-');
hold on;
box on;
grid on;
semilogy(X2,Y2,'DisplayName',"$\beta_R = 5.0$",'LineWidth',5,'Marker','o','MarkerSize',20,'LineStyle','-');
semilogy(X3,Y3,'DisplayName',"$\beta_R = 1.0$",'LineWidth',5,'Marker','x','MarkerSize',20,'LineStyle','-');
semilogy(X4,Y4,'DisplayName',"$\beta_R = 0.5$",'LineWidth',5,'Marker','+','MarkerSize',20,'LineStyle','-');
xlabel("$s$",'Interpreter','latex');
ylabel("$\tau^{\rm pop,coh}_{\rm opt}$",'Interpreter','latex')
fontsize(gca,45,"pixels")
legend('location','southeast','Interpreter','latex')
yline(1e-6,'--','Linewidth',5,'HandleVisibility','off');

xlim([0.75 3]);

%set(gca,'Ytick',logspace(-5,-1,7));
set(gca, 'YScale', 'log', 'YTickLabelMode', 'auto');
hold off;

savefig(sprintf('coh_plot_new2_NL1=%d,e=%.2f,g=%.4f,ham_type=%d.fig',NL1,e,g,ham_type))