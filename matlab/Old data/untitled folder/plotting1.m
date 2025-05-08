%% plotting nicely... NL+1

%X1=getfield(load('./data_1_1.mat'),'betalist');
Y1=getfield(load('./coh_data_1.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
Y2=getfield(load('./coh_data_2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
Y3=getfield(load('./coh_data_3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
Y4=getfield(load('./coh_data_4.mat'),'cvx_optval');

Y5=getfield(load('./coh_data_5.mat'),'cvx_optval');
Y6=getfield(load('./coh_data_6.mat'),'cvx_optval');
Y7=getfield(load('./coh_data_7.mat'),'cvx_optval');
Y8=getfield(load('./coh_data_8.mat'),'cvx_optval');
Y9=getfield(load('./coh_data_9.mat'),'cvx_optval');
Y10=getfield(load('./coh_data_10.mat'),'cvx_optval');

%X1=getfield(load('./data_1_1.mat'),'betalist');
Z1=getfield(load('./coh_data2_1.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
Z2=getfield(load('./coh_data2_2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
Z3=getfield(load('./coh_data2_3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
Z4=getfield(load('./coh_data2_4.mat'),'cvx_optval');

Z5=getfield(load('./coh_data2_5.mat'),'cvx_optval');
Z6=getfield(load('./coh_data2_6.mat'),'cvx_optval');
Z7=getfield(load('./coh_data2_7.mat'),'cvx_optval');
Z8=getfield(load('./coh_data2_8.mat'),'cvx_optval');
Z9=getfield(load('./coh_data2_9.mat'),'cvx_optval');
Z10=getfield(load('./coh_data2_10.mat'),'cvx_optval');

%X1=getfield(load('./data_1_1.mat'),'betalist');
W1=getfield(load('./coh_data3_1.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
W2=getfield(load('./coh_data3_2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
W3=getfield(load('./coh_data3_3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
W4=getfield(load('./coh_data3_4.mat'),'cvx_optval');

W5=getfield(load('./coh_data3_5.mat'),'cvx_optval');
W6=getfield(load('./coh_data3_6.mat'),'cvx_optval');
W7=getfield(load('./coh_data3_7.mat'),'cvx_optval');
W8=getfield(load('./coh_data3_8.mat'),'cvx_optval');
W9=getfield(load('./coh_data3_9.mat'),'cvx_optval');
W10=getfield(load('./coh_data3_10.mat'),'cvx_optval');

%X1=getfield(load('./data_1_1.mat'),'betalist');
V1=getfield(load('./coh_data4_1.mat'),'cvx_optval');

%X2=getfield(load('./data_1_2.mat'),'betalist');
V2=getfield(load('./coh_data4_2.mat'),'cvx_optval');

%X3=getfield(load('./data_1_3.mat'),'betalist');
V3=getfield(load('./coh_data4_3.mat'),'cvx_optval');

%X4=getfield(load('./data_1_4.mat'),'betalist');
V4=getfield(load('./coh_data4_4.mat'),'cvx_optval');

V5=getfield(load('./coh_data4_5.mat'),'cvx_optval');
V6=getfield(load('./coh_data4_6.mat'),'cvx_optval');
V7=getfield(load('./coh_data4_7.mat'),'cvx_optval');
V8=getfield(load('./coh_data4_8.mat'),'cvx_optval');
V9=getfield(load('./coh_data4_9.mat'),'cvx_optval');
V10=getfield(load('./coh_data4_10.mat'),'cvx_optval');
%% make sure all these elists are the same !!! Makes comparison easier!!

X_1 = [0.001,0.00215443,0.00464159,0.01,0.02154435,0.04641589,0.1,0.21544347,0.46415888,1];
Y_1 = [Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10];

Z_1 = [Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10];

W_1 = [W1,W2,W3,W4,W5,W6,W7,W8,W9,W10];
V_1 = [V1,V2,V3,V4,V5,V6,V7,V8,V9,V10];



loglog(X_1,Y_1,'DisplayName',"$N_M = 0$",'LineWidth',3,'Marker','*','MarkerSize',15,'LineStyle','-');
hold on;
box on;
grid on;
loglog(X_1,Z_1,'DisplayName',"$N_M = 1$",'LineWidth',3,'Marker','o','MarkerSize',15,'LineStyle','-');
loglog(X_1,W_1,'DisplayName',"$N_M = 2$",'LineWidth',3,'Marker','x','MarkerSize',15,'LineStyle','-');
loglog(X_1,V_1,'DisplayName',"$N_M = 3$",'LineWidth',3,'Marker','+','MarkerSize',15,'LineStyle','-');
xlabel("$g$",'Interpreter','latex');
ylabel("$\tau_{opt coh}$",'Interpreter','latex')
fontsize(gca,45,"pixels")
legend('location','northwest','Interpreter','latex')
yline(1e-6,'--','Linewidth',3,'HandleVisibility','off');

%set(gca,'Ytick',logspace(-12,-4,8));
hold off;

savefig('plot1.fig')