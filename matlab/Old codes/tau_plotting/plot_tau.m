% plotting nicely... NL = 2

X1=getfield(load('./data_1.mat'),'fixed_glist');
Y1=getfield(load('./data_1.mat'),'tau_list');

X2=getfield(load('./data_2.mat'),'fixed_glist');
Y2=getfield(load('./data_2.mat'),'tau_list');

X3=getfield(load('./data_3.mat'),'fixed_glist');
Y3=getfield(load('./data_3.mat'),'tau_list');

X4=getfield(load('./data_4.mat'),'fixed_glist');
Y4=getfield(load('./data_4.mat'),'tau_list');

X5=getfield(load('./data_5.mat'),'fixed_glist');
Y5=getfield(load('./data_5.mat'),'tau_list');

%% make sure all these elists are the same !!! Makes comparison easier!!



loglog(X1,Y1,'DisplayName',"(i)",'LineWidth',3,'Marker','*','MarkerSize',15,'LineStyle','-');

hold on;
grid on;
box on;

loglog(X2,Y2,'DisplayName',"(ii)",'LineWidth',3,'Marker','o','MarkerSize',15,'LineStyle','-');
loglog(X3,Y3,'DisplayName',"(iii)",'LineWidth',3,'Marker','x','MarkerSize',15,'LineStyle','-');
loglog(X4,Y4,'DisplayName',"(iv)",'LineWidth',3,'Marker','+','MarkerSize',15,'LineStyle','-');
loglog(X5,Y5,'DisplayName',"(v)",'LineWidth',3,'Marker','s','MarkerSize',15,'LineStyle','-');
xlabel("$g_4$",'Interpreter','latex');
ylabel("$\tau$",'Interpreter','latex');
fontsize(gca,45,"pixels")

legend('location','eastoutside','fontsize',30)

savefig('plot3.fig')