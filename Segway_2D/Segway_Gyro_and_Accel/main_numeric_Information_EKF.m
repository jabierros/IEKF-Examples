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
% Segway Gyro and Accel
R=0.25
delta1=0.0
delta3=1
m_R=15
m_P=80
I_R=1/2*(0.5*m_R)*R^2
I_P=m_P*(delta1^2+delta3^2)/(4)^2
c=(m_R+m_P)*0.75/2 %Delta_s=v*t+1/2*a*t^2=2m=2m/s*4s+1/2*a*(4s)^2 => a=(2-8)/(8)=0.75
g=9.8

Delta_t=0.001

param=[R delta1 delta3 m_R m_P I_R I_P c g Delta_t ]' % as is from main_symbolic_EKF.m
%%---Edit_Begin

% Set Initial time
%%---Edit_Begin
t_0=0;       %s
%%---Edit_End

% Actual initial state x_actual_0
%%---Edit_Begin
x_0=0
theta_0=0.01*pi/180
dx_0=0
dtheta_0=0

q_0=[x_0,theta_0]'
dq_0=[dx_0,dtheta_0]'
M_m=0
x_actual_0=[q_0;dq_0]
%%---Edit_End

%Initial mu_x covariance (assumed diagonal)
%%---Edit_Begin
sigma_x0=10 % m
sigma_theta0=30*pi/180 % rad
sigma_dx0=1 %  m/s
sigma_dtheta0=1/(delta1^2+delta3^2)^0.5 %angular velocity if at  1m/s the wheel stops abruptly
diag_Sigma_x_0=[sigma_x0,sigma_theta0,sigma_dx0,sigma_dtheta0]';
%%---Edit_End

% Initial mu_x
%%---Edit_Begin
mu_x_0 = x_actual_0+diag_Sigma_x_0; % you can leave it as is
%%---Edit_End

% Parameters for input measurement noise
%%---Edit_Begin
sigma_M_m_actual=0.1; %Nm
sigma_M_m_spec=sigma_M_m_actual; %Nm
Sigma_u_actual=diag([sigma_M_m_actual])
diag_Sigma_u=[sigma_M_m_spec]; %input meas. cov. assumed diagonal
%%---Edit_End

% Parameters for process discretization noise
%%---Edit_Begin
max_error_discr=[1.8278e-5 , 1.8278e-5 ,0.0004295870 ,0.0004295870]'
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
sigma_gyro_actual=0.5; %rad/s
sigma_accx_actual=0.1*10; %m/s^2
sigma_accz_actual=0.1*10; %m/s^2
%Sensor error and sensor error siven in spec. sheet are equal.
sigma_gyro_spec=sigma_gyro_actual; %rad/s
sigma_accx_spec=sigma_accx_actual; %m/s^2
sigma_accz_spec=sigma_accz_actual; %m/s^2

Sigma_z_actual=diag([sigma_gyro_actual,sigma_accx_actual,sigma_accz_actual]);
diag_Sigma_z=[sigma_gyro_spec,sigma_accx_spec,sigma_accz_spec]; %sensor measurement error assumed diagonal
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
q_string=["x";
"\theta"];
dq_string=["\dot{x}";
"\dot{\theta}"];
x_string=[q_string;
          dq_string];

u_string=["M^m"];

z_string=["\omega_y";
          "a_x";
          "a_z"]
%---End_Edit

%% Kalman Filter Loop
%Starting simulation time and length
%%---Edit_Begin
t_0=0;
t_end=15
%%---Edit_End

%% Maximum Likelihood theta_=[diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u] determination
use_previously_found_ML_params=true
if use_previously_found_ML_params
%%---Edit_Begin
    theta_ =[1.8461e-6    2.7149e-6    0.2807e-6  361.2978e-6    0.4969    1.0050    0.9927    0.1417 ]';
%%---Edit_End
    diag_Sigma_discr=theta_(1:n_x)
    diag_Sigma_z=theta_(n_x+1:n_x+n_z)
    diag_Sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u)
    logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0)
end

ML_find_params=false
if ML_find_params
theta_=[diag_Sigma_discr;diag_Sigma_z;diag_Sigma_u];
% function parameters: logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0)
fun = @(theta_) logL_IEKF(theta_(1:n_x),theta_(n_x+1:n_x+n_z),theta_(n_x+n_z+1:n_x+n_z+n_u),diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0);
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_,options);

diag_Sigma_discr=theta_(1:n_x)
diag_Sigma_z=theta_(n_x+1:n_x+n_z)
diag_Sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u)
end

%% Maximum Likelihood theta_=[diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,mu_x_0,diag_Sigma_x_0] determination
use_previously_found_ML_params=false;
if use_previously_found_ML_params
%%---Edit_Begin
    theta_ =[2.6624e-6   10.9100e-6    1.0400e-6  400.6525e-6    0.4972    1.0053    0.9905   -0.2474    8.4982    0.0013    0.4821    0.1302   16.7719    0.0013    0.5878    0.1294]';
%%---Edit_End
    diag_Sigma_discr=theta_(1:n_x)
    diag_Sigma_z=theta_(n_x+1:n_x+n_z)
    diag_Sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u)
    mu_x_0=theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x)
    diag_Sigma_x_0=theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x)
    logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0)
end

%%
ML_find_params=false
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

num_samples_statistic=1000;

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
% M_m=0
u_actual=25*(heaviside(t-5)-heaviside(t-10));
end
%%---Edit_End
