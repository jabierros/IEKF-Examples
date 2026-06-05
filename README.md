# IEKF Examples

This repository contains a generalized **MATLAB** implementation of the Information Extended Kalman Filter (IEKF), Information Extended Kalman Smoother (IEKS), and Iterated IEKS (IIEKS). It includes Maximum Likelihood Estimation (MLE) for optimizing filter parameters.

---

## 1. Project Directory Structure

* **`LibIEKF/`**: Contains the core filter and smoother pass managers:
  * `IEKF.m` / `IEKF_step.m`: The forward pass and individual epoch filter updates.
  * `IEKS.m` / `IEKS_step.m`: The backward pass and individual smoother updates.
  * `IIEKS.m`: Iterated smoothing.
  * `Plotting.m` / `unpack_simulation.m`: Log analysis and visualization.
  * `logL_IEKF.m`: Negative log-likelihood and joint energy computation.
* **`LibIEKF/Template/`**: Generic templates for the simulator interface:
  * `get_u.m` / `get_z.m` / `get_x_true.m`: Integrate continuous equations using `ode45` and generate sensor/input measurements contaminated with noise.
* **`Double_Pendulum/`**: A 2-Degree-of-Freedom (MDOF) double pendulum example.
* **`Mass_Spring_Damper/`**: A 1-Degree-of-Freedom (1D) mass-spring-damper example.

---

## 2. Decoupled Struct-Based Architecture

Rather than relying on global variables, `LibIEKF` separates system parameters, simulation options, and physical state representations into three structures:

1. **`KF` (Filter Settings)**:
   * `Delta_t`: Filter time-step.
   * `param`: Model parameters vector passed to analytical functions.
   * `mu_x_0` / `sigma_x_0`: Initial filter state mean and standard deviation.
   * `Sigma_w` / `Sigma_v` / `Sigma_u` / `Sigma_z`: Cholesky noise covariance factors.
2. **`TrueSystem` (Truth Simulation)**:
   * `param_true`: True physical system parameters.
   * `x_true_0`: Starting physical state vector.
   * `sigma_u_true` / `sigma_z_true`: Noise levels of inputs and sensors.
3. **`SimOpts` (Simulation Options)**:
   * `t_0` / `t_end`: Start and end times.
   * `t` / `t_prev`: Current and previous time epochs.

---

## 3. How to Run

Before running the examples, ensure that the library and template folders are added to your MATLAB path:
```matlab
addpath('LibIEKF');
addpath('LibIEKF/Template');
```

For either example folder:
1. **Regenerate Equations**: Run `main_symbolic_EKF.m` inside the example directory. This derives and exports analytical transition, measurement, and Jacobian functions ending with trailing underscores (`f_`, `h_`, `f_x_`, etc.) to prevent namespace conflicts.
2. **Run Simulation & Filter**: Run `main_numeric_Information_EKF.m` to simulate the true trajectory, execute the filter and smoother passes, and perform parameter identification via MLE.

---

## 4. Typical Invocation & Log-Likelihood Optimization

### Filter & Smoother Passes
```matlab
% Staged datalogging keys
datalogging_string={'t';'x_true';'mu_x';'Sigma2_x';'u_true';'u_meas';'z_true';'z_meas';'sigma_x'};

% Run forward IEKF pass
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);

% Run backward IEKS smoother pass
FilterResults = IEKS(FilterResults, KF, SimOpts);
```

### Parameter Tuning via MLE
To tune filter noise parameters (like process noise std `sigma_w`), compute the prediction negative log-likelihood or smoother joint energy and optimize using `fminsearch`:
```matlab
% Set up initial noise parameters
theta = [sigma_w];

% Define anonymous objective function
fun = @(theta) logL_IEKF(update_struct(KF, 'sigma_w', theta), TrueSystem, SimOpts);

% Run optimization
options = optimset('PlotFcns', @optimplotfval);
theta = fminsearch(fun, theta, options);
```

---
*The content below corresponds to the older version (`LibIEKF_0`) of the library, which relied on global variables:*

## [Old Version] Globals
Extensive use of global variables is made on purpose, just to keep the code as simple as possible.

* `global t t_0 t_end Delta_t`. time, initial time, final time and time step length (they are same for filter and simulation).
* `global mu_x Sigma2_x`. state vector and state covariance matrix determined by the filter.
* `global u_meas z_meas`. simulated measurement vector of filter input and sensors. `u_meas=get_u()` and `z_meas=get_z()`, are called to get the variables at time `t`. When calling these functions the actual system state is integrated to actualize it to the current time.
* `global sigma_u_actual sigma_z_actual x_actual x_actual_0 u_actual_func seed`.  Actual system measurement covariance of input and sensor (make `sigma_u_actual=zero(n_u,n_u)` if input isn't noisy), actual system state, actual system initial state. `u_actual_func` is a function handle `u_actual_func = @(t) (...)` to a function defining the actual system input as a function of time. Other global variables in this epigraph (actual system variables) can be used to implement the function. `seed` actual system random generator seed used to generate measurements of input and sensors. As the input is generated alongside the filter, parameter tuning requires repeatable random measurement sequences, and to that end seed is used to restart the random number generator at each (`t_0:Delta_t:t_end`) invocation of the filter.
* `global param` these are the model and sensor equation parameters (coincident with those defined in `main_symbolic_EKF.m`)

## [Old Version] Typical invocation
### Filter alone
```
t_0=...;
t=t_0;
Delta_t=...;
t_end=...;% t_end>t_0+Delta_t

mu_x=...; n_x=size(x,1);
Sigma2_x =...; if not(isequal(size(Sigma2_x),[n_x,n_x])) || not(issymmetric(Sigma2_x)) || not(all(eig(Sigma2_x) >  0) error('Sigma2_x is n_x times n_x, symmetric' and positive definite); end

Sigma_u_actual=..;
Sigma_z_actual=..;
x_actual=x_actual_0;
seed=1789;
u_actual_func =...; % u_actual_func = @(t) my_u_actual_func(t);

sigma_discr=...; n_x times 1 column vector of filter assumed discretization error std-s
sigma_z=...; n_z times 1 column vector of filter assumed sensor measurement std-s
sigma_u=...; n_u times 1 column vector of filter assumed input measurement std-s
sigma_w_x=...; n_x times 1 column vector of filter assumed other process equation noise (related to modeling errors usually)
sigma_v_x=...; n_z times 1 column vector of filter assumed other sensor equation noise (related to modeling errors usually)

rng(seed); u_meas=get_u(); z_meas=get_z();

for k=1:t_end/Delta_t
    IEKF(sigma_discr, sigma_z, sigma_u,sigma_w_x,sigma_v_x);
end
```

### Filter iteration with datalogging
```
...
fid=fopen('sol.dat','w'); % Sets file for datalogging

datalogging_string={'t';'x_actual';'mu_x';'u_actual';'u_meas';'z_actual';'z_meas';'sigma_x';'mu_x_error'}; % Set variables for dataloging

% Datalog k=0 %Requires restarting seed
datalogging(fid, datalogging_string);

rng(seed); u_meas=get_u(); z_meas=get_z();

for k=1:t_end/Delta_t
    % filter
    IEKF(sigma_discr, sigma_z, sigma_u,sigma_w_x,sigma_v_x);

    % Datalog k+1
    datalogging(fid, datalogging_string)

end

fclose(fid);
```

### Retrieving log and plotting

```
load_datalogging('sol.dat', datalogging_string);
t_series, x_actual_series, mu_x_series, u_actual_series, u_meas_series, z_actual_series, z_meas_series, sigma_x_series
mu_x_error_series=mu_x_series-x_actual_series;
plot(t_series,mu_x_series-x_actual_series,'-');
```

### Determine Linear Observability
```
load_datalogging('sol.dat', datalogging_string)
k=size(t_series,1) % linear observability near the end of the simulation
f_x_=f_x(mu_x_series(k),u_meas_series(k),t_series(k),param);
h_x_=h_x(mu_x_series(k),u_meas_series(k),t_series(k),param);
clear *_series

OB=obsv(f_x_,h_x_)
rank(OB)
```

### Determine a kind of Nonlinear Observability
```
load_datalogging('sol.dat', datalogging_string)
k=size(t_series,1) % nonlinear observability near the end of the simulation
k=k-5

OB=[h_x(x_actual_series(k,:)',u_meas_series(k,:)',t_series(k,:),param)
    h_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)*f_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)
    h_x(x_actual_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_actual_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)
    h_x(x_actual_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x(x_actual_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x(x_actual_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)
    h_x(x_actual_series(k+4,:)',u_meas_series(k+4,:)',t_series(k+4,:),param)*f_x(x_actual_series(k+4,:)',u_meas_series(k+4,:)',t_series(k+4,:),param)*f_x(x_actual_series(k+3,:)',u_meas_series(k+3,:)',t_series(k+3,:),param)*f_x(x_actual_series(k+2,:)',u_meas_series(k+2,:)',t_series(k+2,:),param)*f_x(x_actual_series(k+1,:)',u_meas_series(k+1,:)',t_series(k+1,:),param)]
clear *_series

rank(OB)
```

### Determine -LogLikelihood of the filter prediction series.
```
logL_IEKF(sigma_discr,sigma_z,sigma_u,sigma_w_x,sigma_v_x,mu_x_0,sigma_x_0)
```

### Minimize -LogLikelihood to identify filter parameters in `theta_` (excluded filter initial state)

```
theta_=[sigma_discr;sigma_z;sigma_u;mu_x_0;sigma_x_0]; % define parameters to be identified
fun = @(theta_) logL_IEKF(theta_(1:2),theta_(3),theta_(4),sigma_w_x,sigma_v_x,theta_(5:6),theta_(7:8));
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_,options);

% untangle parameters
sigma_discr=theta_(1:2)
sigma_z=theta_(3)
sigma_u=theta_(4)
mu_x_0=theta_(5:6)
sigma_x_0=theta_(7:8)
```

### Minimize -LogLikelihood to identify filter parameters `theta_` (including filter initial state)
```
theta_=[sigma_discr;sigma_z;sigma_u;mu_x_0;sigma_x_0]; % set the vector of to-be-identified parameters
fun = @(theta_) logL_IEKF(theta_(1:2),theta_(3),theta_(4),sigma_w_x,sigma_v_x,theta_(5:6),theta_(7:8));
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_,options);

% untangle parameters
sigma_discr=theta_(1:2)
sigma_z=theta_(3)
sigma_u=theta_(4)
mu_x_0=theta_(5:6)
sigma_x_0=theta_(7:8)
```
