clear all; clc;
% cd(fileparts(which(mfilename)));  % In a .m file uncomment this
cd(fileparts(matlab.desktop.editor.getActiveFilename)); % In a .mlx file uncomment this
%% Example: Double Pendulum

%% Solution:
% Parameters and initial state definition

%% 
% Define handle to function determining the actual input

u_true_func = @u_true_func;
%% 

seed=1789; rng(seed);
rng_status=rng;
%% 

% Set Initial time

t_0=0;         %s
Delta_t=0.001; %s
t_end=25;      %s
%% 
% Define model parameters

m1=1;         %kg
m2=1;         %kg
l1=1;         %m
l2=1;         %m
g=10;         %m/s^2
param=[m1 m2 l1 l2 g Delta_t]'; % as is from main_symbolic_EKF.m
param_true=param;
%% 

theta1_0=pi/6;      %rad
dtheta1_0=0;        %rad/s
theta2_0=pi/6;      %rad
dtheta2_0=0;        %rad/s
q_0=[theta1_0,theta2_0]';
dq_0=[dtheta1_0,dtheta2_0]';
x_true_0=[q_0;dq_0];
%% 

sigma_theta1_0=pi/2;  % rad no idea at all
sigma_dtheta1_0=0.1;  % rad/s no idea at all
sigma_theta2_0=pi/2;  % rad no idea at all
sigma_dtheta2_0=0.1;  % rad/s no idea at all
sigma_x_0= [sigma_theta1_0,sigma_dtheta1_0,sigma_theta2_0,sigma_dtheta2_0]';

mu_x_0 = x_true_0+sigma_x_0;
Sigma2_x_0 = diag(sigma_x_0.^2); % assumed diagonal @ t=t_0
%% 
% Model equation error variance $\sigma^2_{u_{k}}$ and $\sigma^2_{w_k}$

n_u=size(u_true_func(t_0),1);
sigma_u_true=zeros(n_u,1);
sigma_u=zeros(n_u,1); %input meas. cov. assumed diagonal

%% 

% max_error_discr=[1/2*1,1/2*1,1/2*1,1/2*1]'*Delta_t^2;
max_error_discr=[2.58418e-06, 3.025e-06, 1.08183e-05, 2.12142e-05]';
%% 
% Parameters for other noise sources in process equation $\mathbf{w}_{k}$.  
% Process equation error std $\mathbf{\sigma}_{\mathbf{w}_{k}}$.

n_x=size(x_true_0,1);
sigma_w_x=max_error_discr; % Discretization error is the process noise source in this example
%% 

sigma_gyroy_true=0.1;%rad/s
sigma_gyroy_spec=sigma_gyroy_true;
sigma_accx_true=0.1;%m/s^2
sigma_accx_spec=sigma_accx_true;
sigma_accz_true=0.1;%m/s^2
sigma_accz_spec=sigma_accz_true;

sigma_z_true=[sigma_gyroy_true,sigma_accx_true,sigma_accz_true]';
sigma_z=[sigma_gyroy_spec,sigma_accx_spec,sigma_accz_spec]';
%% 
% Parameters for other noise sources in sensor equation $\mathbf{v}_{k+1}$. 
% Sensor equation error std $\mathbf{\sigma}_{\mathbf{v}_{k+1}}$.

n_z=size(sigma_z_true,1);
sigma_v_x=zeros(n_z,1);

% Define and initialize the state struct S
S = struct();
S.t_0 = t_0;
S.Delta_t = Delta_t;
S.t_end = t_end;
S.t = t_0;
S.t_prev = t_0;
S.param = param;
S.param_true = param_true;
S.x_true_0 = x_true_0;
S.x_true = x_true_0;
S.x_true_prev = x_true_0;
S.mu_x_0 = mu_x_0;
S.mu_x = mu_x_0;
S.sigma_x_0 = sigma_x_0;
S.Sigma2_x_0 = Sigma2_x_0;
S.Sigma2_x = Sigma2_x_0;
S.sigma_u_true = sigma_u_true;
S.sigma_z_true = sigma_z_true;
S.sigma_u = sigma_u;
S.sigma_z = sigma_z;
S.sigma_w_x = sigma_w_x;
S.sigma_v_x = sigma_v_x;
S.u_true_func = u_true_func;

datalogging_string={'t';'x_true';'mu_x';'Sigma2_x';'u_true';'u_meas';'z_true';'z_meas';'sigma_x'};

%% Kalman Filter Loop (Initial)
% Reset random number generator and run simulation
rng(rng_status);
Simulation = IEKF(S, datalogging_string);
Simulation = IEKS(Simulation, S);
% Simulation = IIEKS(Simulation, S, 4);

%% Observability
% Linear observability near the end of the simulation
[u_meas_tmp, S] = get_u(S); [z_meas_tmp, S] = get_z(S);
f_x_=f_x(S.mu_x, u_meas_tmp, S.t, S.param);
h_x_=h_x(S.mu_x, u_meas_tmp, S.t, S.param);
 
OB=obsv(f_x_,h_x_);
rank(OB)
size(OB) 
%% 
% Nonlinear observability near the end of the simulation
unpack_simulation(Simulation);
%% 
% Nonlinear observabilty matrix must be bigger (more rows) than the linear one. 
% So 5 time steps must suffice ( $5 \times 2=10>8$ )

k=length(t_series)-5;

OB=[h_x(x_true_series(k,:)',u_meas_series(k,:)',t_series(k,:),param)
    h_x(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)*f_x(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)
    h_x(x_true_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_true_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)
    h_x(x_true_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x(x_true_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x(x_true_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)
    h_x(x_true_series(k+4,:)',u_meas_series(k+4,:)',t_series(k+4,:),param)*f_x(x_true_series(k+4,:)',u_meas_series(k+4,:)',t_series(k+4,:),param)*f_x(x_true_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x(x_true_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)];
rank(OB)
clear *_series
%% Real filter error statistics in the limit when $k \longrightarrow \infty$

unpack_simulation(Simulation);
mu_x_error_series=mu_x_series-x_true_series;
num_samples_statistic=100;
lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;
%% Plotting

q_string=["\theta_1";"\theta_2"];
dq_string=["\dot{\theta}_1";"\dot{\theta}_2"];
x_string=[q_string; dq_string];
u_string=[];
z_string=["\omega_y";"a_x";"a_z"];

S.x_string = x_string;
S.u_string = u_string;
S.z_string = z_string;

unpack_simulation(Simulation);
fig_dir=['IEKF_',num2str(Delta_t)];
Plotting(fig_dir, Simulation, S);

%% Stop the code so user takes control
return;

%% Maximum Likelihood (ML) Estimation of Filter Parameters (Opt 1) - Run Optimization
% Run this section manually to search for parameters.

theta_=[S.sigma_w_x];
fun = @(theta_) logL_IEKF(update_struct(S, 'sigma_w_x', theta_));
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_, options);
title('$LogML$','Interpreter','latex');
xlabel('Iteration','Interpreter','latex');
ylabel('$LogML$ value','Interpreter','latex');

sigma_w_x=theta_(1:n_x);
sigma_z=S.sigma_z;
sigma_u=S.sigma_u;
save('ML_opt_1','sigma_w_x','sigma_z','sigma_u');

%% Load Optimization Results & Run Filter (Opt 1)
load('ML_opt_1');
S.sigma_w_x = sigma_w_x;
S.sigma_z = sigma_z;
S.sigma_u = sigma_u;

%% Kalman Filter Loop (Opt 1)
rng(rng_status);
Simulation = IEKF(S, datalogging_string);
Simulation = IEKS(Simulation, S);
% Simulation = IIEKS(Simulation, S, 4);

%% Real filter error statistics (Opt 1)
unpack_simulation(Simulation);
mu_x_error_series=mu_x_series-x_true_series;
num_samples_statistic=100;
lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

%% Plotting (Opt 1)
Plotting(fig_dir, Simulation, S);

%% Stop the code so user takes control
return;

%% Maximum Likelihood (ML) Estimation of Filter Parameters (Opt 2) - Run Optimization
% Run this section manually to search for parameters.

theta_=[S.sigma_w_x; S.sigma_z; S.sigma_u];
fun = @(theta_) logL_IEKF(update_struct(S, ...
    'sigma_w_x', theta_(1:n_x), ...
    'sigma_z', theta_(n_x+1:n_x+n_z), ...
    'sigma_u', theta_(n_x+n_z+1:n_x+n_z+n_u)));
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_, options);

sigma_w_x=theta_(1:n_x);
sigma_z=theta_(n_x+1:n_x+n_z);
sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u);

title('$LogML$','Interpreter','latex');
xlabel('Iteration','Interpreter','latex');
ylabel('$LogML$ value','Interpreter','latex');
save('ML_opt_2','sigma_w_x','sigma_z','sigma_u');

%% Load Optimization Results & Run Filter (Opt 2)
load('ML_opt_2');
S.sigma_w_x = sigma_w_x;
S.sigma_z = sigma_z;
S.sigma_u = sigma_u;

%% Kalman Filter Loop (Opt 2)
rng(rng_status);
Simulation = IEKF(S, datalogging_string);
Simulation = IEKS(Simulation, S);
% Simulation = IIEKS(Simulation, S, 4);

%% Real filter error statistics (Opt 2)
unpack_simulation(Simulation);
mu_x_error_series=mu_x_series-x_true_series;
num_samples_statistic=100;
lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

%% Plotting (Opt 2)
Plotting(fig_dir, Simulation, S);

%% Stop the code so user takes control
return;

%% Maximum Likelihood (ML) Estimation of Filter Parameters (Opt 3) - Run Optimization
% Run this section manually to search for parameters.

theta_=[S.sigma_w_x; S.sigma_z; S.sigma_u; S.mu_x_0; S.sigma_x_0];
fun = @(theta_) logL_IEKF(update_struct(S, ...
    'sigma_w_x', theta_(1:n_x), ...
    'sigma_z', theta_(n_x+1:n_x+n_z), ...
    'sigma_u', theta_(n_x+n_z+1:n_x+n_z+n_u), ...
    'mu_x_0', theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x), ...
    'sigma_x_0', theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x)));
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_, options);
title('$LogML$','Interpreter','latex');
xlabel('Iteration','Interpreter','latex');
ylabel('$LogML$ value','Interpreter','latex');

sigma_w_x=theta_(1:n_x);
sigma_z=theta_(n_x+1:n_x+n_z);
sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u);
mu_x_0=theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x);
sigma_x_0=theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x);
save('ML_opt_3','sigma_w_x','sigma_z','sigma_u','mu_x_0','sigma_x_0');

%% Load Optimization Results & Run Filter (Opt 3)
load('ML_opt_3');
S.sigma_w_x = sigma_w_x;
S.sigma_z = sigma_z;
S.sigma_u = sigma_u;
S.mu_x_0 = mu_x_0;
S.sigma_x_0 = sigma_x_0;
S.Sigma2_x_0 = diag(sigma_x_0.^2);

%% Kalman Filter Loop (Opt 3)
rng(rng_status);
Simulation = IEKF(S, datalogging_string);
Simulation = IEKS(Simulation, S);
% Simulation = IIEKS(Simulation, S, 4);

%% Real filter error statistics (Opt 3)
unpack_simulation(Simulation);
mu_x_error_series=mu_x_series-x_true_series;
num_samples_statistic=100;
lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

%% Plotting (Opt 3)
Plotting(fig_dir, Simulation, S);

%% Stop the code so user takes control
return;

%% Maximum Likelihood (ML) Estimation of Filter & Model Parameters (Opt 4) - Run Optimization
% Run this section manually to search for parameters.

theta_=[S.sigma_w_x; S.sigma_z; S.sigma_u; S.mu_x_0; S.sigma_x_0; S.param];
fun = @(theta_) logL_IEKF(update_struct(S, ...
    'sigma_w_x', theta_(1:n_x), ...
    'sigma_z', theta_(n_x+1:n_x+n_z), ...
    'sigma_u', theta_(n_x+n_z+1:n_x+n_z+n_u), ...
    'mu_x_0', theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x), ...
    'sigma_x_0', theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x), ...
    'param', theta_(n_x+n_z+n_u+n_x+n_x+1 : end)));
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_, options);

n_param = length(S.param);
sigma_w_x=theta_(1:n_x);
sigma_z=theta_(n_x+1:n_x+n_z);
sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u);
mu_x_0=theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x);
sigma_x_0=theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x);
param=theta_(n_x+n_z+n_u+n_x+n_x+1 : n_x+n_z+n_u+n_x+n_x+n_param);

title('$LogML$ with Model Params','Interpreter','latex');
xlabel('Iteration','Interpreter','latex');
ylabel('$LogML$ value','Interpreter','latex');
save('ML_opt_4','sigma_w_x','sigma_z','sigma_u','mu_x_0','sigma_x_0','param');

%% Load Optimization Results & Run Filter (Opt 4)
load('ML_opt_4');
S.sigma_w_x = sigma_w_x;
S.sigma_z = sigma_z;
S.sigma_u = sigma_u;
S.mu_x_0 = mu_x_0;
S.sigma_x_0 = sigma_x_0;
S.Sigma2_x_0 = diag(sigma_x_0.^2);
S.param = param;

%% Kalman Filter Loop (Opt 4)
rng(rng_status);
Simulation = IEKF(S, datalogging_string);
Simulation = IEKS(Simulation, S);
% Simulation = IIEKS(Simulation, S, 4);

%% Real filter error statistics (Opt 4)
unpack_simulation(Simulation);
mu_x_error_series=mu_x_series-x_true_series;
num_samples_statistic=100;
lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

%% Plotting (Opt 4)
Plotting(fig_dir, Simulation, S);




function S = update_struct(S, varargin)
    for i = 1:2:length(varargin)
        S.(varargin{i}) = varargin{i+1};
    end
end