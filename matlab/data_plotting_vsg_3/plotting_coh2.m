e = 0.01;
ham_type=2;
NL1=2;

X1=getfield(load('./coh_data_NL1=3,NM=0.mat'),'fixed_glist');
Y1=getfield(load('./coh_data_NL1=3,NM=0.mat'),'tau_list');

X2=getfield(load('./coh_data_NL1=3,NM=1.mat'),'fixed_glist');
Y2=getfield(load('./coh_data_NL1=3,NM=1.mat'),'tau_list');

X3=getfield(load('./coh_data_NL1=3,NM=2.mat'),'fixed_glist');
Y3=getfield(load('./coh_data_NL1=3,NM=2.mat'),'tau_list');

X4=getfield(load('./coh_data_NL1=3,NM=3.mat'),'fixed_glist'); 
Y4=getfield(load('./coh_data_NL1=3,NM=3.mat'),'tau_list');

%% make sure all these elists are the same !!! Makes comparison easier!!




loglog(X1,Y1,'DisplayName',"$N_M = 0$",'LineWidth',3,'Marker','*','MarkerSize',15,'LineStyle','-');
hold on;
box on;
grid on;
loglog(X2,Y2,'DisplayName',"$N_M = 1$",'LineWidth',3,'Marker','o','MarkerSize',15,'LineStyle','-');
loglog(X3,Y3,'DisplayName',"$N_M = 2$",'LineWidth',3,'Marker','x','MarkerSize',15,'LineStyle','-');
loglog(X4,Y4,'DisplayName',"$N_M = 3$",'LineWidth',3,'Marker','+','MarkerSize',15,'LineStyle','-');
xlabel("$g$",'Interpreter','latex');
ylabel("$\tau^{coh}_{opt}$",'Interpreter','latex')
fontsize(gca,45,"pixels")
legend('location','northwest','Interpreter','latex')
yline(1e-6,'--','Linewidth',3,'HandleVisibility','off');

set(gca,'Ytick',logspace(-14,-2,13));
hold off;

savefig('coh_plot_NL1=3.fig')