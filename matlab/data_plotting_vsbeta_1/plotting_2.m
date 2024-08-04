% plotting nicely... NL = 2

X1=getfield(load('./data_2_1.mat'),'betalist');
Y1=getfield(load('./data_2_1.mat'),'optimal_value');

X2=getfield(load('./data_2_2.mat'),'betalist');
Y2=getfield(load('./data_2_2.mat'),'optimal_value');

X3=getfield(load('./data_2_3.mat'),'betalist');
Y3=getfield(load('./data_2_3.mat'),'optimal_value');

X4=getfield(load('./data_2_4.mat'),'betalist');
Y4=getfield(load('./data_2_4.mat'),'optimal_value');

X5=getfield(load('./data_2_5.mat'),'betalist');
Y5=getfield(load('./data_2_5.mat'),'optimal_value');

%% make sure all these elists are the same !!! Makes comparison easier!!






loglog(X1,Y1,'DisplayName',"$N_M = 1$",'LineWidth',3,'Marker','*','MarkerSize',15,'LineStyle','-');
hold on;
box on;
grid on;
loglog(X2,Y2,'DisplayName',"$N_M = 2$",'LineWidth',3,'Marker','o','MarkerSize',15,'LineStyle','-');
loglog(X3,Y3,'DisplayName',"$N_M = 3$",'LineWidth',3,'Marker','x','MarkerSize',15,'LineStyle','-');
loglog(X4,Y4,'DisplayName',"$N_M = 4$",'LineWidth',3,'Marker','+','MarkerSize',15,'LineStyle','-');
loglog(X5,Y5,'DisplayName',"$N_M = 5$",'LineWidth',3,'Marker','s','MarkerSize',15,'LineStyle','-');

xlabel("$\beta$",'Interpreter','latex');
ylabel("$\tau_{opt}$",'Interpreter','latex')
fontsize(gca,45,"pixels")
%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
legend('location','northwest','Interpreter','latex')
yline(1e-6,'--','Linewidth',3,'HandleVisibility','off');

%set(gca,'Ytick',logspace(-5,-2,4));
savefig('plot2.fig')
