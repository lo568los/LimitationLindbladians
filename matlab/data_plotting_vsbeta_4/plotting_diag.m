e = 0.01;
ham_type=1;
g = 0.0016;
NL1=1;


X1=getfield(load(sprintf('./diag_data_NL1=%d,e=%.2f,beta_r=10.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'betal_list');
Y1=getfield(load(sprintf('./diag_data_NL1=%d,e=%.2f,beta_r=10.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'optimal_value');

X2=getfield(load(sprintf('./diag_data_NL1=%d,e=%.2f,beta_r=5.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'betal_list');
Y2=getfield(load(sprintf('./diag_data_NL1=%d,e=%.2f,beta_r=5.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'optimal_value');

X3=getfield(load(sprintf('./diag_data_NL1=%d,e=%.2f,beta_r=1.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'betal_list');
Y3=getfield(load(sprintf('./diag_data_NL1=%d,e=%.2f,beta_r=1.0,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'optimal_value');

X4=getfield(load(sprintf('./diag_data_NL1=%d,e=%.2f,beta_r=0.5,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'betal_list');
Y4=getfield(load(sprintf('./diag_data_NL1=%d,e=%.2f,beta_r=0.5,g=%.4f,ham_type=%d.mat',NL1,e,g,ham_type)),'optimal_value');

%% make sure all these elists are the same !!! Makes comparison easier!!




loglog(X1,Y1,'DisplayName',"$\beta_R = 10.0$",'LineWidth',5,'Marker','*','MarkerSize',20,'LineStyle','-');
hold on;
box on;
grid on;
loglog(X2,Y2,'DisplayName',"$\beta_R = 5.0$",'LineWidth',5,'Marker','o','MarkerSize',20,'LineStyle','-');
loglog(X3,Y3,'DisplayName',"$\beta_R = 1.0$",'LineWidth',5,'Marker','x','MarkerSize',20,'LineStyle','-');
loglog(X4,Y4,'DisplayName',"$\beta_R = 0.5$",'LineWidth',5,'Marker','+','MarkerSize',20,'LineStyle','-');
xlabel("$\beta_L$",'Interpreter','latex');
ylabel("$\tau_{\rm opt}$",'Interpreter','latex')
fontsize(gca,45,"pixels")
legend('location','southwest','Interpreter','latex')
yline(1e-6,'--','Linewidth',5,'HandleVisibility','off');

set(gca,'Ytick',logspace(-8,0,9));
hold off;

savefig(sprintf('diag_plot_NL1=%d,e=%.2f,g=%.4f,ham_type=%d.fig',NL1,e,g,ham_type))