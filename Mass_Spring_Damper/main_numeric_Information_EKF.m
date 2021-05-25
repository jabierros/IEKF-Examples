clear all
close all

global t t_0 t_end Delta_t
global mu_x Sigma2_x n_x
global u_meas z_meas n_u n_z          
global Sigma_z_actual Sigma_u_actual x_actual x_actual_0 u_actual_func seed % for 'reality' simulation
global param % model and sensor equation parameters (defined in main_symbolic_EKF.m)

% Define actual handle to function determining the input
u_actual_func = @(t) my_u_actual_func(t);

%% Parameters and initial state definition
%%---Edit_Begin
% Mass Spring Damper with ddx sensor
% Set Model and sensor param vector
m=1          %kg
k=20         %Kg/s^2
rho0=1         %m
c=0.1        %Kg/s

Delta_t=0.01 %s

param=[m k rho0 c Delta_t]' % as is from main_symbolic_EKF.m
%%---Edit_End

% Set Initial time
%%---Edit_Begin
t_0=0;       %s
%%---Edit_End

% Actual initial state x_actual_0
%%---Edit_Begin
x_0=2        %m
dx_0=0       %m/s

q_0=[x_0]'
dq_0=[dx_0]'

x_actual_0=[q_0;dq_0]
%%---Edit_End

%Initial mu_x covariance (assumed diagonal)
%%---Edit_Begin
sigma_x0=0.5  % m no idea at all
sigma_dx0=sigma_x0*sqrt(k/m) % m/s no idea at all
diag_Sigma_x_0= [sigma_x0,sigma_dx0]';
%%---Edit_End

% Initial mu_x
%%---Edit_Begin
mu_x_0 = x_actual_0+diag_Sigma_x_0; % you can leave it as is
%%---Edit_End

% Parameters for input measurement noise
%%---Edit_Begin
sigma_f_ext_actual=0.01;   %actual std of the sensor meassuring f_ext, bias assumed 0
sigma_f_ext_spec=sigma_f_ext_actual; %actual std of the sensor meassuring f_ext, bias assumed 0
Sigma_u_actual=diag([sigma_f_ext_actual])
diag_Sigma_u=[sigma_f_ext_spec]; %input meas. cov. assumed diagonal
%%---Edit_End

% Parameters for process discretization noise
%%---Edit_Begin
max_error_discr=[0.001,0.00439535];
diag_Sigma_discr=max_error_discr; %discr. err. cov. assumed diagonal
%%---Edit_End

% Parameters for other noise sources in process quation w_x
%%---Edit_Begin
sigma_w_x_1=0;
diag_Sigma_w_x=[sigma_w_x_1]; %other proccess noise cov. assumed diagonal
%%---Edit_End

% Parameters for other noise sources in process quation v_x
%%---Edit_Begin
sigma_v_x_1=0;
diag_Sigma_v_x=[sigma_v_x_1]; %other proccess  noise cov. assumed diagonal
%%---Edit_End

% Parameters for sensor measurement noise
%%---Edit_Begin
sigma_accx_actual=0.1;
%Sensor error and sensor error siven in spec. sheet are equal.
sigma_accx_spec=sigma_accx_actual;

Sigma_z_actual=diag([sigma_accx_actual]);
diag_Sigma_z=[sigma_accx_spec]; %sensor measurement error assumed diagonal
%%---Edit_End

% Choose random number generator seed
%%---Edit_Begin
seed=1789;
%%---Edit_End

t=t_0;
x_actual=x_actual_0
mu_x=mu_x_0;
Sigma2_x = diag(diag_Sigma_x_0.^2);

n_x=size(x_actual_0,1);
n_u=size(Sigma_u_actual,1);
n_z=size(Sigma_z_actual,1);

%% Kalman Filter Loop
%Starting simulation time and length
%%---Edit_Begin
t_0=0;
t_end=50;
%%---Edit_End

%% Maximum Likelihood theta_=[diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u] determination
use_previously_found_ML_params=false
if use_previously_found_ML_params
%%---Edit_Begin
    theta_ =[19.957691456003325e-4, 7.240187142211141e-4, 0.093495231579440, 0.005939902115621];
%%---Edit_End
    diag_Sigma_discr=theta_(1:n_x)
    diag_Sigma_z=theta_(n_x+1:n_x+n_z)
    diag_Sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u)
    logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0)
end

ML_find_params=false
if ML_find_params
theta_=[diag_Sigma_discr;diag_Sigma_z;diag_Sigma_u];
% function parameters logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,mu_x_0,diag_Sigma_x_0)
fun = @(theta_) logL_IEKF(theta_(1:2),theta_(3),theta_(4),mu_x_0,diag_Sigma_x_0);
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_,options);

diag_Sigma_discr=theta_(1:n_x)
diag_Sigma_z=theta_(n_x+1:n_x+n_z)
diag_Sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u)
end

%% Maximum Likelihood theta_=[diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,mu_x_0,diag_Sigma_x_0] determination
use_previously_found_ML_params=true;
if use_previously_found_ML_params
%%---Edit_Begin
    theta_ =1.0e+02 *[   0.000001830224861e-4, 1.334642513215480e-4, 0.000913294982475, 0.000380941571911, 0.020000027551784, -0.000836123483162, 0.000000054530327, 0.000000348279808]'
%%---Edit_End
    diag_Sigma_discr=theta_(1:n_x)
    diag_Sigma_z=theta_(n_x+1:n_x+n_z)
    diag_Sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u)
    mu_x_0=theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x)
    diag_Sigma_x_0=theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x)
    logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0)
end

%%
ML_find_params=true
if ML_find_params
theta_=[diag_Sigma_discr;diag_Sigma_z;diag_Sigma_u;mu_x_0;diag_Sigma_x_0];
% function parameters: logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0)
fun = @(theta_) logL_IEKF(theta_(1:n_x),theta_(n_x+1:n_x+n_z),theta_(n_x+n_z+1:n_x+n_z+n_u),diag_Sigma_w_x,diag_Sigma_v_x,theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x),theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x));
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_,options);

    diag_Sigma_discr=theta_(1:n_x)
    diag_Sigma_z=theta_(n_x+1:n_x+n_z)
    diag_Sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u)
    mu_x_0=theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x)
    diag_Sigma_x_0=theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x)
end

%% Run filter for t_0:t_end while Datalog Datalogging
% Read Input and Sensor at first (k=1)

% Init filter state
t=t_0;
mu_x=mu_x_0;
Sigma2_x=diag((diag_Sigma_x_0).^2);
rng(seed); u_meas=get_u(); z_meas=get_z();

% Datalogging Setup
% Set file for datalogging
fid=fopen('sol.dat','w');

% Set variables for dataloging
datalogging_string={'t';'x_actual';'mu_x';'u_actual';'u_meas';'z_actual';'z_meas';'diag_Sigma_x'};

% Datalog k=0 %Requires restarting seed
datalogging(fid, datalogging_string);

for k=1:t_end/Delta_t
    % filter
    IEKF(diag_Sigma_discr, diag_Sigma_z, diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x);
    
    % Datalog k+1
    datalogging(fid, datalogging_string)
end

fclose(fid);

%% Observability
% Linear observability near the end of the simulation
f_x_=f_x(mu_x,u_meas,t,param);
h_x_=h_x(mu_x,u_meas,t,param);
 
OB=obsv(f_x_,h_x_)
rank(OB)
size(OB) %Nonlinear observabilty matrix must be bigger (more rows) than the linear one.

% Nonlinear observability near the end of the simulation

load_datalogging('sol.dat', datalogging_string)

size(OB) %Nonlinear observabilty matrix must be bigger (more rows) than the linear one.
k=k-5

OB=[h_x(x_actual_series(k,:)',u_meas_series(k,:)',t_series(k,:),param)
    h_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)*f_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)
    h_x(x_actual_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_actual_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)
    h_x(x_actual_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x(x_actual_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x(x_actual_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)
    h_x(x_actual_series(k+4,:)',u_meas_series(k+4,:)',t_series(k+4,:),param)*f_x(x_actual_series(k+4,:)',u_meas_series(k+4,:)',t_series(k+4,:),param)*f_x(x_actual_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x(x_actual_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)]

rank(OB)

clear *_series

%% Real filter error statistics @ lim k -> infty
load_datalogging('sol.dat', datalogging_string)

mu_x_error_series=mu_x_series-x_actual_series;

num_samples_statistic=100;

lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error)
lim_mu_x_error_std=std(lim_mu_x_error)
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5


%% Plotting
load_datalogging('sol.dat', datalogging_string)

if ML_find_params || use_previously_found_ML_params
    fig_dir=['Information_EKF_ML_find_params_',num2str(Delta_t)];
else
    fig_dir=['Information_EKF_',num2str(Delta_t)];
end
close all
Plotting(fig_dir, datalogging_string)

%% Helper functions
%%
function u_actual=my_u_actual_func(t)
u_actual=0.*t;
end

%%
function logL=logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0)
global t t_0 t_end Delta_t
global u_meas z_meas 
global mu_x Sigma2_x
global mu_x_pred Sigma2_x_pred R
global seed
global param

%% Kalman Filter Loop

% Starting simulation time and length
logL=0;
%% Read Input and Sensor at first (k=1)

rng(seed);

t=t_0;
mu_x=mu_x_0;
Sigma2_x=diag((diag_Sigma_x_0).^2);
rng(seed); u_meas=get_u(); z_meas=get_z();

for k=1:t_end/Delta_t

    IEKF(diag_Sigma_discr, diag_Sigma_z, diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x);
    
    %% Log-Likelihood computation
    h_x_=h_x(mu_x_pred,u_meas,t,param);
    Sigma2_Innovation=h_x_*Sigma2_x_pred*h_x_'+R;
    Innovation=z_meas-h(mu_x_pred,u_meas,t,param);
    logL=logL+1/2*(log(det(Sigma2_Innovation))+Innovation'*inv(Sigma2_Innovation)*Innovation);
end
end

function IEKF(diag_Sigma_discr, diag_Sigma_z, diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x)
global t  Delta_t
global u_meas z_meas
global mu_x Sigma2_x 
global mu_x_pred Sigma2_x_pred Q R
global param
    
    %---Begin--- One step of Kalman filter (k-th)
    
    % Prediction
    u_meas=get_u();
    % Effect of noise in u on model covariance (f_u_ is dependent on the state)
    f_u_=f_u(mu_x,u_meas,t,param);
    Q_u= f_u_*diag(diag_Sigma_u.^2)*f_u_';
    
    % Effect of noise source w_x_ on model covariance (f_w_x_ is dependent on the state)
    f_w_x_=f_w_x(mu_x,u_meas,t,param);
    Q_w_x= f_w_x_*diag(diag_Sigma_w_x.^2)*f_w_x_';
    
    % Discretization error effect on model covariance
    Q_discr=diag((diag_Sigma_discr).^2);
    
    % Error sources are assumed independent
    Q = Q_discr+Q_u+Q_w_x;
    
    f_x_=f_x(mu_x,u_meas,t,param);
    
    mu_x_pred=f(mu_x,u_meas,t,param);
    Sigma2_x_pred=f_x_*Sigma2_x*f_x_'+Q;
    I_pred=inv(Sigma2_x_pred);
    i_pred=I_pred*mu_x_pred;
    
    % Read Input and Sensor (k+1)
    t_prev=t;
    t=t+Delta_t;
    
    z_meas=get_z();
    
    % Projection / Observation
    h_u_= h_u(mu_x_pred,u_meas,t,param);
    % effect of noise in u on sensor equation covariance
    R_u= h_u_*diag(diag_Sigma_u.^2)*h_u_';
    
    h_v_x_= h_v_x(mu_x_pred,u_meas,t,param) ;
    % effect of noise source w_x_ on sensor equation covariance
    R_v_x= h_v_x_*diag(diag_Sigma_v_x.^2)*h_v_x_';
    
    R_z=diag(diag_Sigma_z.^2);
    % Sensor error, input error, and other noise sources (w_x_,...) error are assumed independent
    R=R_z+R_u+R_v_x;
    
    h_x_=h_x(mu_x_pred,u_meas,t,param);
    i_proj=h_x_'*inv(R)*(z_meas-(h(mu_x_pred,u_meas,t,param)-h_x_*mu_x_pred));
    I_proj=h_x_'*inv(R)*h_x_;
    
    % Fussion / Information Weighted Average
    I=I_pred+I_proj;
    %mu_x=inv(I)*(I_pred*mu_x_pred+I_proj*mu_x_proj);
    mu_x=inv(I)*(i_pred+i_proj);
    Sigma2_x=inv(I);
    
    %%---End--- One step of Kalman filter (k-th)

end

%%
function set_x_actual()
global t t_0
global x_actual_0
global x_actual u_actual_func
global param
persistent t_prev x_actual_prev
if t==t_0
    x_actual=x_actual_0;
    x_actual_prev=x_actual_0;
    t_prev=t_0;
elseif t==t_prev
    ;
elseif t>t_prev
    %True system state is generated by integration using MATLAB ode45 integrator
    ode45_options=odeset('RelTol',1e-12,'AbsTol',1e-12);
    [time_steps,x_steps] = ode45(@(t,x_) dstate(x_, u_actual_func(t_prev), t_prev,param),[t_prev,t],x_actual_prev,ode45_options);
    x_actual=x_steps(end,:)';
    x_actual_prev=x_actual;
    t_prev=t;
else
    error('t must be such t==t_0 or t=t_prev+Delta_t');
end
end

function u_meas=get_u()
global t
global u_actual
global Sigma_u_actual u_actual_func 
set_x_actual()

u_actual=u_actual_func(t);
n_u=size(u_actual,1);
u_meas=u_actual+Sigma_u_actual*randn(n_u,1);
end

function z_meas=get_z()
global t
global x_actual z_actual u_actual u_actual_func 
global Sigma_z_actual
global param
set_x_actual()
u_actual=u_actual_func(t);

%Sensor data is generated from true system estate and input (x,u), by using
%h, to which sensor noise is added
z_actual=h(x_actual,u_actual,t,param);
n_z=size(z_actual,1);
z_meas=z_actual+Sigma_z_actual*randn(n_z,1);

end

%%
function datalogging(fid, datalogging_string)
% Datalogged variables
global t u_meas u_actual z_meas z_actual x_actual mu_x Sigma2_x diag_Sigma_x

diag_Sigma_x=diag(Sigma2_x).^0.5;

for i=1:size(datalogging_string,1)
    var=eval(datalogging_string{i});
    for j=1:size(var,1)
     fprintf(fid,'%d ',var(j,1));
    end
end
fprintf(fid,'\n');
end

%%
function load_datalogging(sol, datalogging_string)
% Datalogged variables
global t u_meas u_actual z_meas z_actual x_actual mu_x Sigma2_x diag_Sigma_x 

load(sol,'-ascii')

column_end=0;
for i=1:size(datalogging_string,1)
    n_columns=size(eval(datalogging_string{i}),1);
    column_start=column_end+1;
    column_end=column_end+n_columns;
    assignin('caller',[datalogging_string{i},'_series'], sol(:,column_start:column_end));
end
clear sol
end

%% Plotting
function Plotting(fig_dir, datalogging_string)
global t_end
global x_actual_0
global param

[mkdir_success,mkdir_message]=mkdir(fig_dir);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

fontsize=12

load_datalogging('sol.dat', datalogging_string);

mu_x_error_series=mu_x_series-x_actual_series;

figIdx=1;

fig=figure(figIdx);
figname='x_actual';
plot(t_series,x_actual_series);
set(gca, 'YScale', 'linear');
legH=legend('${x}$','$\dot{x}$');
titH=title('$\mathbf{x}$ (Actual State)');
set(fig,'units','normalized');
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize)
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

x_discr_series=[];
x_=x_actual_0;
for k=1:length(t_series)
 x_discr_series=[x_discr_series,x_];
 x_=f(x_, u_actual_series(k,:)', t_series(k),param);
end

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='x_discr';
plot(t_series,x_discr_series);
set(gca, 'YScale', 'linear');
legH=legend('${x}$','$\dot{x}$');
titH=title('$\mathbf{x}$ (integrated with discr.)');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='u_meas_u_actual';
plot(t_series,u_meas_series,'-',t_series,u_actual_series,'--');
set(gca, 'YScale', 'linear');
legH=legend('$\mathrm{meas}(f^{ext})$', '$f^{ext}$');
titH=title('$\mathrm{meas}(\mathbf{u})$ vs $\mathbf{u}$ ');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='u_meas_minus_u_actual';
plot(t_series,u_meas_series-u_actual_series,'-');
set(gca, 'YScale', 'linear');
legH=legend('$\mathrm{meas}(f^{ext})-f^{ext}$')
titH=title('$\mathrm{meas}(\mathbf{z})-\mathbf{z}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='z_meas_z_actual';
plot(t_series,z_meas_series,'-',t_series,z_actual_series,'--');
set(gca, 'YScale', 'linear');
legH=legend('$\mathrm{meas}(\ddot{x})$','$\ddot{x}$');
titH=title('$\mathrm{meas}(\mathbf{z})$ vs $\mathbf{z}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='z_meas_minus_z_actual';
plot(t_series,z_meas_series-z_actual_series,'-');
set(gca, 'YScale', 'linear');
legH=legend('$\mathrm{meas}(\ddot{x})-\ddot{x}$');
titH=title('$\mathrm{meas}(\mathbf{z})-\mathbf{z}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='mu_x_x_actual';
plot(t_series,mu_x_series,'-',t_series,x_actual_series,'--');
set(gca, 'YScale', 'linear');
xlim([0 t_end]);
ylim([-6 6]);
legH= legend('$\hat{\mu}_x$','$\hat{\mu}_{\dot{x}}$','$x$','$\dot{x}$');
titH=title('$\mathbf{x}$ vs $\hat{{\mu}}_{\mathbf{x}}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

figIdx = figIdx + 1;
fig=figure(figIdx);
set(fig,'units','normalized'); 
plot(t_series,mu_x_series-x_actual_series,'-');
set(gca, 'YScale', 'linear');
xlim([0 t_end]);
ylim([-1 1]);
legH= legend('$\hat{\mu}_{x}-{x}$','$\hat{\mu}_{\dot{x}}-\dot{x}$');
set(legH,'interpreter','latex','Fontsize',fontsize);
titH=title('$\hat{\mu}_{\mathbf{x}}-\mathbf{x}$');
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
figname='mu_x_minus_x_actual';
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

% Real filter error statistics @ lim k -> infty
num_samples_statistic=100;

lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);;
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='sigma_x';
plot(t_series,diag_Sigma_x_series,'-',t_series,ones(size(diag_Sigma_x_series,1),1)*sqrt_lim_mu_x_error_squared_mean,'--');
set(gca, 'YScale', 'log');
yl = ylim;
ylim([1e-3,yl(2)]);
legH=legend('$\sigma_{x}$','$\sigma_{\dot{x}}$','${\mathrm{mean}((\hat{\mu}_{x}-{x})^2)^{\frac{1}{2}}}$','${\mathrm{mean}((\hat{\mu}_{\dot{x}}-\dot{x})^2)^{\frac{1}{2}}}$');
titH=title('$\mathrm{diag}({{\Sigma}}^2_{\mathbf{x}})^{\frac{1}{2}}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

pred_error_series=zeros(size(x_actual_series,2),size(t_series,1));
for k=1:size(t_series,1)-1
    pred_error_series(:,k+1)=abs((f(x_actual_series(k,:)',u_actual_series(k,:)',t_series(k,:)',param)-x_actual_series(k,:)')-(x_actual_series(k+1,:)'-x_actual_series(k,:)'));
end
pred_error_series=pred_error_series';

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='discr_error';
h=plot(t_series,pred_error_series);
legH=legend('$\varepsilon_{x}$','$\varepsilon_{\dot{x}}$');
titH=title('Discretization error');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

dcm_obj = datacursormode(gcf);
for i=1:length(h)
set(gca, 'YScale', 'log');
[y_max,x_max]=max(pred_error_series(:,i));hDatatip(i)=createDatatip(dcm_obj, h(i), [t_series(x_max),pred_error_series(x_max,i)]);
set(hDatatip(i),'Position',[t_series(x_max),pred_error_series(x_max,i)]);
set(hDatatip(i),'interpreter','latex');
end

set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'epsc');
system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

end


