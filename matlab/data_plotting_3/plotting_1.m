%% plotting nicely... NL+1

X1=getfield(load('./data_1_1.mat'),'fixed_glist');
Y1=getfield(load('./data_1_1.mat'),'optimal_value');

X2=getfield(load('./data_1_2.mat'),'fixed_glist');
Y2=getfield(load('./data_1_2.mat'),'optimal_value');

X3=getfield(load('./data_1_3.mat'),'fixed_glist');
Y3=getfield(load('./data_1_3.mat'),'optimal_value');

X4=getfield(load('./data_1_4.mat'),'fixed_glist');
Y4=getfield(load('./data_1_4.mat'),'optimal_value');
%% make sure all these elists are the same !!! Makes comparison easier!!




loglog(X1,Y1,'DisplayName',"$N_M = 1$",'LineWidth',3,'Marker','*','MarkerSize',15,'LineStyle','-');
hold on;
box on;
grid on;
loglog(X2,Y2,'DisplayName',"$N_M = 2$",'LineWidth',3,'Marker','o','MarkerSize',15,'LineStyle','-');
loglog(X3,Y3,'DisplayName',"$N_M = 3$",'LineWidth',3,'Marker','x','MarkerSize',15,'LineStyle','-');
loglog(X4,Y4,'DisplayName',"$N_M = 4$",'LineWidth',3,'Marker','+','MarkerSize',15,'LineStyle','-');
xlabel("$g$",'Interpreter','latex');
ylabel("$\tau_{opt}$",'Interpreter','latex')
fontsize(gca,45,"pixels")
legend('location','northwest','Interpreter','latex')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
hold off;

savefig('plot1.fig')