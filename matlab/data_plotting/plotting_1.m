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


hold on;
box on;
grid on;
plot(X1,Y1,'DisplayName',"N_M = 1",'LineWidth',3,'Marker','*','MarkerSize',15,'LineStyle','-');
plot(X2,Y2,'DisplayName',"N_M = 2",'LineWidth',3,'Marker','o','MarkerSize',15,'LineStyle','-');
plot(X3,Y3,'DisplayName',"N_M = 3",'LineWidth',3,'Marker','x','MarkerSize',15,'LineStyle','-');
plot(X4,Y4,'DisplayName',"N_M = 3",'LineWidth',3,'Marker','+','MarkerSize',15,'LineStyle','-');
xlabel("g");
ylabel("\tau_{opt}")
fontsize(gca,42,"pixels")
legend()

hold off;

savefig('plot1.fig')