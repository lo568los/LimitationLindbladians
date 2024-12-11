

%X1=getfield(load('./data_1_1.mat'),'betalist');
W1=getfield(load('./thermal_data_1.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
W2=getfield(load('./thermal_data_2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
W3=getfield(load('./thermal_data_3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
W4=getfield(load('./thermal_data_4.mat'),'cvx_optval');

W5=getfield(load('./thermal_data_5.mat'),'cvx_optval');
W6=getfield(load('./thermal_data_6.mat'),'cvx_optval');
W7=getfield(load('./thermal_data_7.mat'),'cvx_optval');
W8=getfield(load('./thermal_data_8.mat'),'cvx_optval');
W9=getfield(load('./thermal_data_9.mat'),'cvx_optval');
W10=getfield(load('./thermal_data_10.mat'),'cvx_optval');

%X1=getfield(load('./data_1_1.mat'),'betalist');
V1=getfield(load('./thermal_data2_1.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
V2=getfield(load('./thermal_data2_2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
V3=getfield(load('./thermal_data2_3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
V4=getfield(load('./thermal_data2_4.mat'),'cvx_optval');

V5=getfield(load('./thermal_data2_5.mat'),'cvx_optval');
V6=getfield(load('./thermal_data2_6.mat'),'cvx_optval');
V7=getfield(load('./thermal_data2_7.mat'),'cvx_optval');
V8=getfield(load('./thermal_data2_8.mat'),'cvx_optval');
V9=getfield(load('./thermal_data2_9.mat'),'cvx_optval');
V10=getfield(load('./thermal_data2_10.mat'),'cvx_optval');

X_1 = [1.0,1.2,1.5,1.7,2.0,2.2,2.5,2.7,3.0,3.5];

W_1 = [W1,W2,W3,W4,W5,W6,W7,W8,W9,W10];
V_1 = [V1,V2,V3,V4,V5,V6,V7,V8,V9,V10];

loglog(X_1,W_1,'DisplayName',"XXZ, e = 0.01",'LineWidth',3,'Marker','x','MarkerSize',15,'LineStyle','-');
hold on;
box on;
grid on;
loglog(X_1,V_1,'DisplayName',"XX, e = 0.01",'LineWidth',3,'Marker','+','MarkerSize',15,'LineStyle','-');
xlabel("$\beta_L$",'Interpreter','latex');
ylabel("$\tau_{opt}$",'Interpreter','latex')
fontsize(gca,45,"pixels")
legend('location','southeast','Interpreter','latex')
yline(1e-6,'--','Linewidth',3,'HandleVisibility','off');

%set(gca,'Ytick',logspace(-12,-4,8));
hold off;

savefig('plot2.fig')