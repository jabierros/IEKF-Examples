clear all
close all

global t t_0 t_end Delta_t
global mu_x Sigma2_x n_x
global u_meas z_meas n_u n_z          
global Sigma_z_actual Sigma_u_actual x_actual x_actual_0 u_actual_func seed % for 'reality' simulation
global param % model and sensor equation parameters (defined in main_symbolic_EKF.m)
global x_string u_string z_string

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

% Variable names for Plotting()
%---Begin_Edit
q_string=["x"];
dq_string=["\dot{x}"];
x_string=[q_string;
          dq_string];

u_string=["f^{ext}"]
z_string=["a_x"];
%---End_Edit

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

Plotting(fig_dir, datalogging_string)

%% Helper functions
%%---Edit_Begin
function u_actual=my_u_actual_func(t)
u_actual=0.*t;
end
%%---Edit_End
