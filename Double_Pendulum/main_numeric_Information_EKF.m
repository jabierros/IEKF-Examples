clear all; clc;
try
    cd(fileparts(matlab.desktop.editor.getActiveFilename));
catch
    try
        cd(fileparts(which(mfilename)));
    catch
        % Keep current directory
    end
end
addpath('../LibIEKF');
addpath('../LibIEKF/Template');
%% Example: Double Pendulum

%% Solution:
% Parameters and initial state definition

seed=1789; rng(seed);
rng_status=rng;

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

n_u=size(u_true_func_(t_0),1);
sigma_u_true=zeros(n_u,1);
sigma_u=zeros(n_u,1); %input meas. cov. assumed diagonal

%% 

% max_error_discr=[1/2*1,1/2*1,1/2*1,1/2*1]'*Delta_t^2;
max_error_discr=[2.58418e-06, 3.025e-06, 1.08183e-05, 2.12142e-05]';
%% 
% Parameters for other noise sources in process equation $\mathbf{w}_{k}$.  
% Process equation error std $\mathbf{\sigma}_{\mathbf{w}_{k}}$.

n_x=size(x_true_0,1);
sigma_w=max_error_discr; % Discretization error is the process noise source in this example
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
sigma_v=zeros(n_z,1);

% Define and initialize the separate configurations and states
SimOpts = struct();
SimOpts.t_0 = t_0;
SimOpts.t_end = t_end;
SimOpts.t = t_0;
SimOpts.t_prev = t_0;

KF = struct();
KF.Delta_t = Delta_t;
KF.param = param;
KF.mu_x_0 = mu_x_0;
KF.sigma_x_0 = sigma_x_0;
KF.Sigma_w = diag(sigma_w);
KF.Sigma_v = diag(sigma_v);
KF.Sigma_u = diag(sigma_u);
KF.Sigma_z = diag(sigma_z);

TrueSystem = struct();
TrueSystem.param_true = param_true;
TrueSystem.x_true_0 = x_true_0;
TrueSystem.x_true = x_true_0;
TrueSystem.x_true_prev = x_true_0;
TrueSystem.sigma_u_true = sigma_u_true;
TrueSystem.sigma_z_true = sigma_z_true;

datalogging_string={'t';'x_true';'mu_x';'Sigma2_x';'u_true';'u_meas';'z_true';'z_meas';'sigma_x'};

%% Kalman Filter Loop (Initial)
% Reset random number generator and run simulation
rng(rng_status);
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);
FilterResults = IEKS(FilterResults, KF, SimOpts);
% FilterResults = IIEKS(FilterResults, KF, SimOpts, 4);

%% Observability
% Linear observability near the end of the simulation
[u_meas_tmp, TrueSystem, SimOpts] = meas_u(TrueSystem, SimOpts); 
[z_meas_tmp, TrueSystem, SimOpts] = meas_z(TrueSystem, SimOpts);
f_x_val=f_x_(FilterResults.data.mu_x{end}, u_meas_tmp, SimOpts.t, KF.param, KF.Delta_t);
h_x_val=h_x_(FilterResults.data.mu_x{end}, u_meas_tmp, SimOpts.t, KF.param);
 
OB=obsv(f_x_val,h_x_val);
rank(OB)
size(OB) 
%% 
% Nonlinear observability near the end of the simulation
unpack_simulation(FilterResults);
%% 
% Nonlinear observabilty matrix must be bigger (more rows) than the linear one. 
% So 5 time steps must suffice ( $5 \times 2=10>8$ )

k=length(t_series)-5;

OB=[h_x_(x_true_series(k,:)',u_meas_series(k,:)',t_series(k,:),param)
    h_x_(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)*f_x_(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param,KF.Delta_t)
    h_x_(x_true_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x_(x_true_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param,KF.Delta_t)*f_x_(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param,KF.Delta_t)
    h_x_(x_true_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x_(x_true_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param,KF.Delta_t)*f_x_(x_true_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param,KF.Delta_t)*f_x_(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param,KF.Delta_t)
    h_x_(x_true_series(k+4,:)',u_meas_series(k+4,:)',t_series(k+4,:),param)*f_x_(x_true_series(k+4,:)',u_meas_series(k+4,:)',t_series(k+4,:),param,KF.Delta_t)*f_x_(x_true_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param,KF.Delta_t)*f_x_(x_true_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param,KF.Delta_t)*f_x_(x_true_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param,KF.Delta_t)];
rank(OB)
clear *_series
%% Real filter error statistics in the limit when $k \longrightarrow \infty$

unpack_simulation(FilterResults);
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

KF.x_string = x_string;
KF.u_string = u_string;
KF.z_string = z_string;

unpack_simulation(FilterResults);
fig_dir=['IEKF_',num2str(Delta_t)];
Plotting(fig_dir, FilterResults, KF, TrueSystem, SimOpts);

%% Stop the code so user takes control
return;

%% Maximum Likelihood (ML) Estimation of Filter Parameters (Opt 1) - Run Optimization
% Run this section manually to search for parameters.

theta_=[sigma_w];
fun = @(theta_) logL_IEKF(update_struct(KF, 'sigma_w', theta_), TrueSystem, SimOpts);
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_, options);
title('$LogML$','Interpreter','latex');
xlabel('Iteration','Interpreter','latex');
ylabel('$LogML$ value','Interpreter','latex');

sigma_w=theta_(1:n_x);
sigma_z=diag(KF.Sigma_z);
sigma_u=diag(KF.Sigma_u);
save('ML_opt_1','sigma_w','sigma_z','sigma_u');

%% Load Optimization Results & Run Filter (Opt 1)
load('ML_opt_1');
if exist('sigma_w_x', 'var') && ~exist('sigma_w', 'var')
    sigma_w = sigma_w_x;
end
KF.Sigma_w = diag(sigma_w);
KF.Sigma_z = diag(sigma_z);
KF.Sigma_u = diag(sigma_u);

%% Kalman Filter Loop (Opt 1)
rng(rng_status);
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);
FilterResults = IEKS(FilterResults, KF, SimOpts);
% FilterResults = IIEKS(FilterResults, KF, SimOpts, 4);

%% Real filter error statistics (Opt 1)
unpack_simulation(FilterResults);
mu_x_error_series=mu_x_series-x_true_series;
num_samples_statistic=100;
lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

%% Plotting (Opt 1)
Plotting(fig_dir, FilterResults, KF, TrueSystem, SimOpts);

%% Stop the code so user takes control
return;

%% Maximum Likelihood (ML) Estimation of Filter Parameters (Opt 2) - Run Optimization
% Run this section manually to search for parameters.

theta_=[sigma_w; diag(KF.Sigma_z); diag(KF.Sigma_u)];
fun = @(theta_) logL_IEKF(update_struct(KF, ...
    'sigma_w', theta_(1:n_x), ...
    'sigma_z', theta_(n_x+1:n_x+n_z), ...
    'sigma_u', theta_(n_x+n_z+1:n_x+n_z+n_u)), TrueSystem, SimOpts);
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_, options);

sigma_w=theta_(1:n_x);
sigma_z=theta_(n_x+1:n_x+n_z);
sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u);

title('$LogML$','Interpreter','latex');
xlabel('Iteration','Interpreter','latex');
ylabel('$LogML$ value','Interpreter','latex');
save('ML_opt_2','sigma_w','sigma_z','sigma_u');

%% Load Optimization Results & Run Filter (Opt 2)
load('ML_opt_2');
if exist('sigma_w_x', 'var') && ~exist('sigma_w', 'var')
    sigma_w = sigma_w_x;
end
KF.Sigma_w = diag(sigma_w);
KF.Sigma_z = diag(sigma_z);
KF.Sigma_u = diag(sigma_u);

%% Kalman Filter Loop (Opt 2)
rng(rng_status);
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);
FilterResults = IEKS(FilterResults, KF, SimOpts);
% FilterResults = IIEKS(FilterResults, KF, SimOpts, 4);

%% Real filter error statistics (Opt 2)
unpack_simulation(FilterResults);
mu_x_error_series=mu_x_series-x_true_series;
num_samples_statistic=100;
lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

%% Plotting (Opt 2)
Plotting(fig_dir, FilterResults, KF, TrueSystem, SimOpts);

%% Stop the code so user takes control
return;

%% Maximum Likelihood (ML) Estimation of Filter Parameters (Opt 3) - Run Optimization
% Run this section manually to search for parameters.

theta_=[sigma_w; diag(KF.Sigma_z); diag(KF.Sigma_u); KF.mu_x_0; KF.sigma_x_0];
fun = @(theta_) logL_IEKF(update_struct(KF, ...
    'sigma_w', theta_(1:n_x), ...
    'sigma_z', theta_(n_x+1:n_x+n_z), ...
    'sigma_u', theta_(n_x+n_z+1:n_x+n_z+n_u), ...
    'mu_x_0', theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x), ...
    'sigma_x_0', theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x)), TrueSystem, SimOpts);
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_, options);
title('$LogML$','Interpreter','latex');
xlabel('Iteration','Interpreter','latex');
ylabel('$LogML$ value','Interpreter','latex');

sigma_w=theta_(1:n_x);
sigma_z=theta_(n_x+1:n_x+n_z);
sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u);
mu_x_0=theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x);
sigma_x_0=theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x);
save('ML_opt_3','sigma_w','sigma_z','sigma_u','mu_x_0','sigma_x_0');

%% Load Optimization Results & Run Filter (Opt 3)
load('ML_opt_3');
if exist('sigma_w_x', 'var') && ~exist('sigma_w', 'var')
    sigma_w = sigma_w_x;
end
KF.Sigma_w = diag(sigma_w);
KF.Sigma_z = diag(sigma_z);
KF.Sigma_u = diag(sigma_u);
KF.mu_x_0 = mu_x_0;
KF.sigma_x_0 = sigma_x_0;

%% Kalman Filter Loop (Opt 3)
rng(rng_status);
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);
FilterResults = IEKS(FilterResults, KF, SimOpts);
% FilterResults = IIEKS(FilterResults, KF, SimOpts, 4);

%% Real filter error statistics (Opt 3)
unpack_simulation(FilterResults);
mu_x_error_series=mu_x_series-x_true_series;
num_samples_statistic=100;
lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

%% Plotting (Opt 3)
Plotting(fig_dir, FilterResults, KF, TrueSystem, SimOpts);

%% Stop the code so user takes control
return;

%% Maximum Likelihood (ML) Estimation of Filter & Model Parameters (Opt 4) - Run Optimization
% Run this section manually to search for parameters.

theta_=[sigma_w; diag(KF.Sigma_z); diag(KF.Sigma_u); KF.mu_x_0; KF.sigma_x_0; KF.param];
fun = @(theta_) logL_IEKF(update_struct(KF, ...
    'sigma_w', theta_(1:n_x), ...
    'sigma_z', theta_(n_x+1:n_x+n_z), ...
    'sigma_u', theta_(n_x+n_z+1:n_x+n_z+n_u), ...
    'mu_x_0', theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x), ...
    'sigma_x_0', theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x), ...
    'param', theta_(n_x+n_z+n_u+n_x+n_x+1 : end)), TrueSystem, SimOpts);
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_, options);

n_param = length(KF.param);
sigma_w=theta_(1:n_x);
sigma_z=theta_(n_x+1:n_x+n_z);
sigma_u=theta_(n_x+n_z+1:n_x+n_z+n_u);
mu_x_0=theta_(n_x+n_z+n_u+1:n_x+n_z+n_u+n_x);
sigma_x_0=theta_(n_x+n_z+n_u+n_x+1:n_x+n_z+n_u+n_x+n_x);
param=theta_(n_x+n_z+n_u+n_x+n_x+1 : n_x+n_z+n_u+n_x+n_x+n_param);

title('$LogML$ with Model Params','Interpreter','latex');
xlabel('Iteration','Interpreter','latex');
ylabel('$LogML$ value','Interpreter','latex');
save('ML_opt_4','sigma_w','sigma_z','sigma_u','mu_x_0','sigma_x_0','param');

%% Load Optimization Results & Run Filter (Opt 4)
load('ML_opt_4');
if exist('sigma_w_x', 'var') && ~exist('sigma_w', 'var')
    sigma_w = sigma_w_x;
end
KF.Sigma_w = diag(sigma_w);
KF.Sigma_z = diag(sigma_z);
KF.Sigma_u = diag(sigma_u);
KF.mu_x_0 = mu_x_0;
KF.sigma_x_0 = sigma_x_0;
KF.param = param;

%% Kalman Filter Loop (Opt 4)
rng(rng_status);
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);
FilterResults = IEKS(FilterResults, KF, SimOpts);
% FilterResults = IIEKS(FilterResults, KF, SimOpts, 4);

%% Real filter error statistics (Opt 4)
unpack_simulation(FilterResults);
mu_x_error_series=mu_x_series-x_true_series;
num_samples_statistic=100;
lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

%% Plotting (Opt 4)
Plotting(fig_dir, FilterResults, KF, TrueSystem, SimOpts);

function KF = update_struct(KF, varargin)
    for i = 1:2:length(varargin)
        name = varargin{i};
        val = varargin{i+1};
        if strcmp(name, 'sigma_w')
            KF.Sigma_w = diag(val);
        elseif strcmp(name, 'sigma_z')
            KF.Sigma_z = diag(val);
        elseif strcmp(name, 'sigma_u')
            KF.Sigma_u = diag(val);
        else
            KF.(name) = val;
        end
    end
end