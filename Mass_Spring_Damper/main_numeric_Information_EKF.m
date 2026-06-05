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

%% Example: Mass Spring Damper

%% Solution:
% Parameters and initial state definition

seed=1789; rng(seed);
rng_status=rng;

% Set initial time $t_0$, sampling frequency $\Delta t$ and final simulation time $t_{end}$
t_0=0;                  % s
Delta_t=0.005;          % s
t_end=50;               % s

% Define model parameters
m=1;          % kg
k=20;         % kg/s^2
rho0=1;       % m
c=0.1;        % kg/s
param=[m k rho0 c]';
param_true=param;

% Set deterministic true/actual system Initial state $x_0^{tr}$
x_0=2;        % m
dx_0=0;       % m/s
q_0=[x_0]';
dq_0=[dx_0]';
x_true_0=[q_0;dq_0];

% Filter Initial state and variance $\mu_{x_0}$, $\sigma_{x_0}^2$, estimated from limits in historical measurements
sigma_x0=0.5;  % m
sigma_dx0=sigma_x0*sqrt(k/m); % m/s
sigma_x_0=[sigma_x0,sigma_dx0]';

mu_x_0 = x_true_0+sigma_x_0;
Sigma2_x_0 = diag(sigma_x_0.^2); % assumed diagonal @ t=t_0

% Model equation error variance $\sigma^2_{u_{k}}$ and $\sigma^2_{w_k}$
n_u=1;
sigma_f_ext=0.01;
sigma_f_ext_spec=0.08;
sigma_u_true=[sigma_f_ext]';
sigma_u=[sigma_f_ext_spec]';

% Discretization error is the process noise source in this example
max_error_discr=[1/2*Delta_t^2,1/2*Delta_t^2];
sigma_w=max_error_discr';

% Sensor equation error std $\sigma_{z_{k+1}}$
sigma_accx_actual=0.1;
sigma_accx_spec=0.8;
sigma_z_true=[sigma_accx_actual]';
sigma_z=[sigma_accx_spec]';

% Other noise sources in sensor equation
n_z=1;
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
[u_meas_tmp, TrueSystem, SimOpts] = get_u(TrueSystem, SimOpts); 
[z_meas_tmp, TrueSystem, SimOpts] = get_z(TrueSystem, SimOpts);
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
q_string=["x"];
dq_string=["\dot{x}"];
x_string=[q_string; dq_string];
u_string=["f^{ext}"];
z_string=["a_x"];

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

sigma_w=theta_(1:size(x_true_0,1));
sigma_z=diag(KF.Sigma_z);
sigma_u=diag(KF.Sigma_u);
save('ML_opt_1','sigma_w','sigma_z','sigma_u');

%% Load Optimization Results & Run Filter (Opt 1)
load('ML_opt_1');
KF.Sigma_w = diag(sigma_w);
KF.Sigma_z = diag(sigma_z);
KF.Sigma_u = diag(sigma_u);

%% Kalman Filter Loop (Opt 1)
rng(rng_status);
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);
FilterResults = IEKS(FilterResults, KF, SimOpts);

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
n_x = size(x_true_0, 1);
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
KF.Sigma_w = diag(sigma_w);
KF.Sigma_z = diag(sigma_z);
KF.Sigma_u = diag(sigma_u);

%% Kalman Filter Loop (Opt 2)
rng(rng_status);
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);
FilterResults = IEKS(FilterResults, KF, SimOpts);

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
n_x = size(x_true_0, 1);
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
KF.Sigma_w = diag(sigma_w);
KF.Sigma_z = diag(sigma_z);
KF.Sigma_u = diag(sigma_u);
KF.mu_x_0 = mu_x_0;
KF.sigma_x_0 = sigma_x_0;

%% Kalman Filter Loop (Opt 3)
rng(rng_status);
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);
FilterResults = IEKS(FilterResults, KF, SimOpts);

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

%% Helper functions
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