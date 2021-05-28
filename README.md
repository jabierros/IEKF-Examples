# IEKF Examples

This is an example of a, pretty general, **MATLAB** implemementation of the Information Extended Kalman Filter.
Maximum Likelihood Estimation of filter parameters is considered.

Real world data is generated algorithmically. To that end an "actual system" state is integrated in parallel with the IEKF. From this "actual" state, "actual" measurements (noise free)  are obtained. These, in turn, are contaminated with noise to generate the so-called "meassured" inputs and sensors. Running in parallel with the IEKF allow to use the filter output as a feedback for control of the system. 

The code is generic enough to be applied to general nonlinear process and sensor equations

## Globals
Extensive use of global variables is made on purpose, just to keep the code as simple as possible.

* `global t t_0 t_end Delta_t`. time, initial time, final time and time step length (they are same for filter and simulation).
* `global mu_x Sigma2_x`. state vector and state covariance matrix determined by the filter.
* `global u_meas z_meas`. simulated measurement vector of filter input and sensors. `u_meas=get_u()` and `z_meas=get_z()`, are called to get the variables at t. When calling these functions the actual system state is integrated to actualize it to the current time.
* `global Sigma_u_actual Sigma_z_actual x_actual x_actual_0 u_actual_func seed`.  Actual system measurement covariance of input and sensor (make `Sigma_u_actual zero(n_u,n_u)` if input isn't noisy), actual system state, initial state. `u_actual_func` is a function handle `u_actual_func = @(t) (...)` to a function defining the actual system input as a function of time. Other global variables in this epigraph (actual sytem variables) can be used to implement the function. `seed` actual system random generator seed used to generate measurements of input and sensors. As the input is generated alongside the filter, parameter tuning requires repeatable random measurement sequences, and to that end seed is used to restart the random number generator at each (`t_0:Delta_t:t_end`) invocation of the filter.
* `global param` these are the model and sensor equation parameters (defined in `main_symbolic_EKF.m`)

## Typical invocation
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

diag_Sigma_discr=...; n_x times 1 column vector of filter assumed discretization error std-s
diag_Sigma_z=...; n_z times 1 column vector of filter assumed sensor measurement std-s
diag_Sigma_u=...; n_u times 1 column vector of filter assumed input measurement std-s
diag_Sigma_w_x=...; n_x times 1 column vector of filter assumed other process equation noise (related to modeling errors usually)
diag_Sigma_v_x=...; n_z times 1 column vector of filter assumed other sensor equation noise (related to modeling errors usually)

rng(seed); u_meas=get_u(); z_meas=get_z();

for k=1:t_end/Delta_t
    IEKF(diag_Sigma_discr, diag_Sigma_z, diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x);
end
```

### Filter iteration with datalogging
```
...
fid=fopen('sol.dat','w'); % Sets file for datalogging

datalogging_string={'t';'x_actual';'mu_x';'u_actual';'u_meas';'z_actual';'z_meas';'diag_Sigma_x';'mu_x_error'}; % Set variables for dataloging

% Datalog k=0 %Requires restarting seed
datalogging(fid, datalogging_string);

rng(seed); u_meas=get_u(); z_meas=get_z();

for k=1:t_end/Delta_t
    % filter
    IEKF(diag_Sigma_discr, diag_Sigma_z, diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x);
    
    % Datalog k+1
    datalogging(fid, datalogging_string)
    
end

fclose(fid);
```

### Retrieving log and plotting

```
load_datalogging('sol.dat', datalogging_string);
t_series, x_actual_series, mu_x_series, u_actual_series, u_meas_series, z_actual_series, z_meas_series, diag_Sigma_x_series
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
logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0)

```

### Minimize -LogLikelihood to identify filter parameters in `theta_` (excluded filter initial state)

```
theta_=[diag_Sigma_discr;diag_Sigma_z;diag_Sigma_u;mu_x_0;diag_Sigma_x_0]; % define parameters to be identified
fun = @(theta_) logL_IEKF(theta_(1:2),theta_(3),theta_(4),diag_Sigma_w_x,diag_Sigma_v_x,theta_(5:6),theta_(7:8));
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_,options);

% untangel parameters
diag_Sigma_discr=theta_(1:2)
diag_Sigma_z=theta_(3)
diag_Sigma_u=theta_(4)
mu_x_0=theta_(5:6)
diag_Sigma_x_0=theta_(7:8)
```

### Minimize -LogLikelihood to identify filter parameters `theta_` (including filter initial state)
```
theta_=[diag_Sigma_discr;diag_Sigma_z;diag_Sigma_u;mu_x_0;diag_Sigma_x_0]; % set the vector of to-be-identified parameters
fun = @(theta_) logL_IEKF(theta_(1:2),theta_(3),theta_(4),diag_Sigma_w_x,diag_Sigma_v_x,theta_(5:6),theta_(7:8));
options = optimset('PlotFcns',@optimplotfval);
theta_ = fminsearch(fun, theta_,options);

% untangle the parameters
diag_Sigma_discr=theta_(1:2)
diag_Sigma_z=theta_(3)
diag_Sigma_u=theta_(4)
mu_x_0=theta_(5:6)
diag_Sigma_x_0=theta_(7:8)

```
