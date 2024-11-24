%% plotting nicely... NL+1

%X1=getfield(load('./data_1_1.mat'),'betalist');
Y1=getfield(load('./thermal_data_new1_th.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
Y2=getfield(load('./thermal_data_new2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
Y3=getfield(load('./thermal_data_new3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
Y4=getfield(load('./thermal_data_new4.mat'),'cvx_optval');

Y5=getfield(load('./thermal_data_new5.mat'),'cvx_optval');
Y6=getfield(load('./thermal_data_new6.mat'),'cvx_optval');
Y7=getfield(load('./thermal_data_new7.mat'),'cvx_optval');
Y8=getfield(load('./thermal_data_new8.mat'),'cvx_optval');
Y9=getfield(load('./thermal_data_new9.mat'),'cvx_optval');
Y10=getfield(load('./thermal_data_new10.mat'),'cvx_optval');

%X1=getfield(load('./data_1_1.mat'),'betalist');
Z1=getfield(load('./thermal_data2_new1.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
Z2=getfield(load('./thermal_data2_new2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
Z3=getfield(load('./thermal_data2_new3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
Z4=getfield(load('./thermal_data2_new4.mat'),'cvx_optval');

Z5=getfield(load('./thermal_data2_new5.mat'),'cvx_optval');
Z6=getfield(load('./thermal_data2_new6.mat'),'cvx_optval');
Z7=getfield(load('./thermal_data2_new7.mat'),'cvx_optval');
Z8=getfield(load('./thermal_data2_new8.mat'),'cvx_optval');
Z9=getfield(load('./thermal_data2_new9.mat'),'cvx_optval');
Z10=getfield(load('./thermal_data2_new10.mat'),'cvx_optval');

%X1=getfield(load('./data_1_1.mat'),'betalist');
W1=getfield(load('./thermal_data3_new1.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
W2=getfield(load('./thermal_data3_new2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
W3=getfield(load('./thermal_data3_new3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
W4=getfield(load('./thermal_data3_new4.mat'),'cvx_optval');

W5=getfield(load('./thermal_data3_new5.mat'),'cvx_optval');
W6=getfield(load('./thermal_data3_new6.mat'),'cvx_optval');
W7=getfield(load('./thermal_data3_new7.mat'),'cvx_optval');
W8=getfield(load('./thermal_data3_new8.mat'),'cvx_optval');
W9=getfield(load('./thermal_data3_new9.mat'),'cvx_optval');
W10=getfield(load('./thermal_data3_new10.mat'),'cvx_optval');

%X1=getfield(load('./data_1_1.mat'),'betalist');
V1=getfield(load('./thermal_data4_new1.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
V2=getfield(load('./thermal_data4_new2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
V3=getfield(load('./thermal_data4_new3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
V4=getfield(load('./thermal_data4_new4.mat'),'cvx_optval');

V5=getfield(load('./thermal_data4_new5.mat'),'cvx_optval');
V6=getfield(load('./thermal_data4_new6.mat'),'cvx_optval');
V7=getfield(load('./thermal_data4_new7.mat'),'cvx_optval');
V8=getfield(load('./thermal_data4_new8.mat'),'cvx_optval');
V9=getfield(load('./thermal_data4_new9.mat'),'cvx_optval');
V10=getfield(load('./thermal_data4_new10.mat'),'cvx_optval');



%% make sure all these elists are the same !!! Makes comparison easier!!

Y_1 = [Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10];
X_1 = [1.0,1.2,1.5,1.7,2.0,2.2,2.5,2.7,3.0,3.5];

Z_1 = [Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10];

W_1 = [W1,W2,W3,W4,W5,W6,W7,W8,W9,W10];
V_1 = [V1,V2,V3,V4,V5,V6,V7,V8,V9,V10];




loglog(X_1,Y_1,'DisplayName',"$\beta_L = 1, N_{L1} = 2$",'LineWidth',3,'Marker','*','MarkerSize',15,'LineStyle','-');
hold on;
box on;
grid on;
loglog(X_1,Z_1,'DisplayName',"$\beta_L = 1.5, N_{L1} = 2$",'LineWidth',3,'Marker','o','MarkerSize',15,'LineStyle','-');
loglog(X_1,W_1,'DisplayName',"$\beta_L = 1, N_{L1} = 1$",'LineWidth',3,'Marker','x','MarkerSize',15,'LineStyle','-');
loglog(X_1,V_1,'DisplayName',"$\beta_L = 1.5, N_{L1} = 1$",'LineWidth',3,'Marker','+','MarkerSize',15,'LineStyle','-');
xlabel("$\beta_R$",'Interpreter','latex');
ylabel("$\tau_{opt}$",'Interpreter','latex')
fontsize(gca,45,"pixels")
legend('location','southeast','Interpreter','latex')
yline(1e-6,'--','Linewidth',3,'HandleVisibility','off');

%set(gca,'Ytick',logspace(-12,-4,8));
hold off;

savefig('plot1.fig')