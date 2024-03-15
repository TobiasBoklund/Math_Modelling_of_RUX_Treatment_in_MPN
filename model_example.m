%% Simulation of model without treatment
% Load parameters
load model_default_parameters.mat;

% Storage of initial conditions
A0 = zeros(8,1);

% Initial conditions
x00 = 1*10^5;
x10 = 2.5*10^6;
x20 = 6.4*10^11;
y00 = 1;
y10 = 0;
y20 = 0;
a0 = 810;
s0 = 1;

% Collect in one vector
A0(1) = x00;
A0(2) = x10;
A0(3) = x20;
A0(4) = y00;
A0(5) = y10;
A0(6) = y20;
A0(7) = a0;
A0(8) = s0;

% Time interval [days]
tmax = 60*365;

% Options
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Solve system of ODEs - use ode15s
[T,A] = ode15s(@model_rhs,[0 tmax],A0,options,p);

% Extract variables
x0 = A(:,1);
x1 = A(:,2);
x2 = A(:,3);
y0 = A(:,4);
y1 = A(:,5);
y2 = A(:,6);
a = A(:,7);
s = A(:,8);

% Calculate allele burden
VAF = y2./(x2+y2);

%% Plots
% Cell plots
fig1 = figure(1);
t1 = tiledlayout(2,2);
xlabel(t1,'$t$/years','fontsize',24,'interpreter','latex');
% fig1.Position = [250 50 1600 1300];
nexttile;
hold on;
grid on;
plot(T/365,x0,'-g','linewidth',6);
plot(T/365,y0,'-r','linewidth',6);
plot(T/365,x0+y0,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',12,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Stem Cells','fontsize',16,'interpreter','latex');
% legend({'$x_0$','$y_0$','$x_0+y_0$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 5;
ylim([0,2*10^5]);
xlim([0,60]);
set(gca,'fontsize',24);

% Legend
leg1 = legend({'Healthy cells','Malignant cells','Sum of healthy and malignant cells'},...
    'orientation','horizontal','interpreter','latex','location','south','fontsize',24);
leg1.Layout.Tile = 'south';

nexttile;
hold on;
grid on;
plot(T/365,x1,'-g','linewidth',6);
plot(T/365,y1,'-r','linewidth',6);
plot(T/365,x1+y1,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Progenitor Cells','fontsize',16,'interpreter','latex');
% legend({'$x_1$','$y_1$','$x_1+y_1$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 6;
ylim([0,8*10^6]);
xlim([0,60]);
set(gca,'fontsize',24);

nexttile;
hold on;
grid on;
plot(T/365,x2,'-g','linewidth',6);
plot(T/365,y2,'-r','linewidth',6);
plot(T/365,x2+y2,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Mature Cells','fontsize',16,'interpreter','latex');
% legend({'$x_2$','$y_2$','$x_2+y_2$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 11;
ylim([10,30*10^11]);
xlim([0,60]);
set(gca,'fontsize',24);

% VAF
nexttile;
hold on;
grid on;
plot(T/365,VAF*100,'-','color','#EDB120','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('\textit{JAK2} VAF/\%','fontsize',14,'interpreter','latex');
title('\textit{JAK2} VAF','fontsize',16,'interpreter','latex');
ylim([0,100]);
xlim([0,60]);
set(gca,'fontsize',24);

% Cell plots
fig2 = figure(2);
t2 = tiledlayout(2,3);
xlabel(t2,'$t$/years','fontsize',24,'interpreter','latex');
% fig1.Position = [250 50 1600 1300];
nexttile;
hold on;
grid on;
plot(T/365,x0,'-g','linewidth',6);
plot(T/365,y0,'-r','linewidth',6);
plot(T/365,x0+y0,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',12,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Stem Cells','fontsize',16,'interpreter','latex');
% legend({'$x_0$','$y_0$','$x_0+y_0$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 5;
ylim([0,2*10^5]);
xlim([0,60]);
set(gca,'fontsize',24);

% Legend
leg1 = legend({'Healthy cells','Malignant cells','Sum of healthy and malignant cells'},...
    'orientation','horizontal','interpreter','latex','location','south','fontsize',24);
leg1.Layout.Tile = 'south';

nexttile;
hold on;
grid on;
plot(T/365,x1,'-g','linewidth',6);
plot(T/365,y1,'-r','linewidth',6);
plot(T/365,x1+y1,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Progenitor Cells','fontsize',16,'interpreter','latex');
% legend({'$x_1$','$y_1$','$x_1+y_1$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 6;
ylim([0,8*10^6]);
xlim([0,60]);
set(gca,'fontsize',24);

nexttile;
hold on;
grid on;
plot(T/365,x2,'-g','linewidth',6);
plot(T/365,y2,'-r','linewidth',6);
plot(T/365,x2+y2,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Mature Cells','fontsize',16,'interpreter','latex');
% legend({'$x_2$','$y_2$','$x_2+y_2$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 11;
ylim([10,30*10^11]);
xlim([0,60]);
set(gca,'fontsize',24);

% VAF
nexttile;
hold on;
grid on;
plot(T/365,VAF*100,'-','color','#EDB120','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('\textit{JAK2} VAF/\%','fontsize',14,'interpreter','latex');
title('\textit{JAK2} VAF','fontsize',16,'interpreter','latex');
ylim([0,100]);
xlim([0,60]);
set(gca,'fontsize',24);

% Cellular debris
nexttile;
hold on;
grid on;
plot(T/365,a,'-','color','#EDB120','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('$a$/1','fontsize',14,'interpreter','latex');
title('Cellular debris','fontsize',16,'interpreter','latex');
ax = gca;
ax.YAxis.Exponent = 3;
ylim([0,2*10^3]);
xlim([0,60]);
set(gca,'fontsize',24);

% Cytokine signal
nexttile;
hold on;
grid on;
plot(T/365,s,'-','color','#EDB120','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('$s$/1','fontsize',14,'interpreter','latex');
title('Cytokine Signal','fontsize',16,'interpreter','latex');
ylim([0,2]);
xlim([0,60]);
set(gca,'fontsize',24);

%% Solve model with treatment
% Load parameters
load model_default_parameters.mat;

% Storage of initial conditions
A0 = zeros(8,1);

% Initial conditions
x00 = 1*10^5;
x10 = 2.5*10^6;
x20 = 6.4*10^11;
y00 = 1;
y10 = 0;
y20 = 0;
a0 = 810;
s0 = 1;

% Collect in one vector
A0(1) = x00;
A0(2) = x10;
A0(3) = x20;
A0(4) = y00;
A0(5) = y10;
A0(6) = y20;
A0(7) = a0;
A0(8) = s0;

% Time interval [days]
tmax = 60*365;

% Define treatment - here set to start after 30 years
ct = linspace(0,tmax,10^5);
cR = zeros(size(ct));
cR(ct>30*365) = 35;

% Define rho values - should be optimised!
rho1 = 0.2;
rho2 = 0.2;

% Define treatment effect
effect = 'sy0dy1';

% Options
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Solve system of ODEs - use ode15s
[T,A] = ode15s(@model_rhs_treatment,[0 tmax],A0,options,p,ct,cR,rho1,rho2,effect);

% Extract variables
x0 = A(:,1);
x1 = A(:,2);
x2 = A(:,3);
y0 = A(:,4);
y1 = A(:,5);
y2 = A(:,6);
a = A(:,7);
s = A(:,8);

% Calculate allele burden
VAF = y2./(x2+y2);

%% Plots
% Cell plots
fig1 = figure(1);
t1 = tiledlayout(2,2);
xlabel(t1,'$t$/years','fontsize',24,'interpreter','latex');
% fig1.Position = [250 50 1600 1300];
nexttile;
hold on;
grid on;
plot(T/365,x0,'-g','linewidth',6);
plot(T/365,y0,'-r','linewidth',6);
plot(T/365,x0+y0,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',12,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Stem Cells','fontsize',16,'interpreter','latex');
% legend({'$x_0$','$y_0$','$x_0+y_0$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 5;
ylim([0,2*10^5]);
xlim([0,60]);
set(gca,'fontsize',24);

% Legend
leg1 = legend({'Healthy cells','Malignant cells','Sum of healthy and malignant cells'},...
    'orientation','horizontal','interpreter','latex','location','south','fontsize',24);
leg1.Layout.Tile = 'south';

nexttile;
hold on;
grid on;
plot(T/365,x1,'-g','linewidth',6);
plot(T/365,y1,'-r','linewidth',6);
plot(T/365,x1+y1,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Progenitor Cells','fontsize',16,'interpreter','latex');
% legend({'$x_1$','$y_1$','$x_1+y_1$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 6;
ylim([0,8*10^6]);
xlim([0,60]);
set(gca,'fontsize',24);

nexttile;
hold on;
grid on;
plot(T/365,x2,'-g','linewidth',6);
plot(T/365,y2,'-r','linewidth',6);
plot(T/365,x2+y2,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Mature Cells','fontsize',16,'interpreter','latex');
% legend({'$x_2$','$y_2$','$x_2+y_2$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 11;
ylim([10,30*10^11]);
xlim([0,60]);
set(gca,'fontsize',24);

% VAF
nexttile;
hold on;
grid on;
plot(T/365,VAF*100,'-','color','#EDB120','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('\textit{JAK2} VAF/\%','fontsize',14,'interpreter','latex');
title('\textit{JAK2} VAF','fontsize',16,'interpreter','latex');
ylim([0,100]);
xlim([0,60]);
set(gca,'fontsize',24);

% Cell plots
fig2 = figure(2);
t2 = tiledlayout(2,3);
xlabel(t2,'$t$/years','fontsize',24,'interpreter','latex');
% fig1.Position = [250 50 1600 1300];
nexttile;
hold on;
grid on;
plot(T/365,x0,'-g','linewidth',6);
plot(T/365,y0,'-r','linewidth',6);
plot(T/365,x0+y0,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',12,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Stem Cells','fontsize',16,'interpreter','latex');
% legend({'$x_0$','$y_0$','$x_0+y_0$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 5;
ylim([0,2*10^5]);
xlim([0,60]);
set(gca,'fontsize',24);

% Legend
leg1 = legend({'Healthy cells','Malignant cells','Sum of healthy and malignant cells'},...
    'orientation','horizontal','interpreter','latex','location','south','fontsize',24);
leg1.Layout.Tile = 'south';

nexttile;
hold on;
grid on;
plot(T/365,x1,'-g','linewidth',6);
plot(T/365,y1,'-r','linewidth',6);
plot(T/365,x1+y1,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Progenitor Cells','fontsize',16,'interpreter','latex');
% legend({'$x_1$','$y_1$','$x_1+y_1$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 6;
ylim([0,8*10^6]);
xlim([0,60]);
set(gca,'fontsize',24);

nexttile;
hold on;
grid on;
plot(T/365,x2,'-g','linewidth',6);
plot(T/365,y2,'-r','linewidth',6);
plot(T/365,x2+y2,'--k','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('Cells/1','fontsize',14,'interpreter','latex');
title('Mature Cells','fontsize',16,'interpreter','latex');
% legend({'$x_2$','$y_2$','$x_2+y_2$'},'fontsize',12,'interpreter','latex',...
%        'location','best');
ax = gca;
ax.YAxis.Exponent = 11;
ylim([10,30*10^11]);
xlim([0,60]);
set(gca,'fontsize',24);

% VAF
nexttile;
hold on;
grid on;
plot(T/365,VAF*100,'-','color','#EDB120','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('\textit{JAK2} VAF/\%','fontsize',14,'interpreter','latex');
title('\textit{JAK2} VAF','fontsize',16,'interpreter','latex');
ylim([0,100]);
xlim([0,60]);
set(gca,'fontsize',24);

% Cellular debris
nexttile;
hold on;
grid on;
plot(T/365,a,'-','color','#EDB120','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('$a$/1','fontsize',14,'interpreter','latex');
title('Cellular debris','fontsize',16,'interpreter','latex');
ax = gca;
ax.YAxis.Exponent = 3;
ylim([0,2*10^3]);
xlim([0,60]);
set(gca,'fontsize',24);

% Cytokine signal
nexttile;
hold on;
grid on;
plot(T/365,s,'-','color','#EDB120','linewidth',6);
% xlabel('$t$/years','fontsize',14,'interpreter','latex');
ylabel('$s$/1','fontsize',14,'interpreter','latex');
title('Cytokine Signal','fontsize',16,'interpreter','latex');
ylim([0,2]);
xlim([0,60]);
set(gca,'fontsize',24);