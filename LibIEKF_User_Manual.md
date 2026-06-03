# LibIEKF User Manual: Writing Symbolic & Numeric Scripts

This manual explains how to write new filter/smoother applications using the **Information Extended Kalman Filter and Smoother (`LibIEKF`)** library. 

---

## 1. Library Architecture & Requirements

`LibIEKF` uses an in-memory, struct-based architecture where all filter states, covariance matrices, parameters, and time steps are managed inside a single configuration structure `S`.

### Required Template Functions
Every system directory (e.g., `Double_Pendulum`) must implement the following three interface functions to supply inputs and measurements:
1. **`get_u.m`**: Fetches or generates the command input `u_meas` (and optionally updates the true system state if doing online simulation).
   ```matlab
   function [u_meas, S] = get_u(S)
   ```
2. **`get_z.m`**: Fetches or generates the measurement vector `z_meas`.
   ```matlab
   function [z_meas, S] = get_z(S)
   ```
3. **`get_x_true.m`**: Integrates continuous physics equations using the true process model (`dstate_true.m`) to generate the true physical trajectory `x_true`.
   ```matlab
   function S = get_x_true(S)
   ```

---

## 2. Writing `main_symbolic.m`

The symbolic script uses MATLAB's Symbolic Math Toolbox to derive the continuous system dynamics, analytical Jacobians, and measurement models, exporting them as optimized MATLAB files.

### Step 1: Define Symbolic Variables & Parameters
Use `syms` to declare coordinates, velocities, time, and physical parameters:
```matlab
syms t theta dtheta ddtheta l m g Delta_t real
x_ = [theta; dtheta];       % Symbolic state vector
u_ = sym(zeros(0,1));       % Symbolic input vector (if none, empty)
param = [m; l; g; Delta_t]; % System parameter vector
```

### Step 2: Formulate System Dynamics
Compute the continuous system derivatives $\mathbf{dstate}$ and discretize using Euler integration:
```matlab
% Continuous state derivative (e.g. from virtual power / Lagrangian dynamics)
dstate = [dtheta; -g/l * sin(theta)]; 

% Discretized process model: x_k+1 = f(x_k, u_k)
f = x_ + dstate * Delta_t;
```

### Step 3: Formulate Measurement Model
Define the measurement function $\mathbf{h}$ in terms of states and inputs:
```matlab
% Example: measuring angular velocity and tangential acceleration
h = [dtheta; -g*sin(theta)]; 
```

### Step 4: Compute Analytical Jacobians
Use `jacobian` to compute the derivative matrices:
```matlab
f_x = jacobian(f, x_);
f_u = jacobian(f, u_);
h_x = jacobian(h, x_);
h_u = jacobian(h, u_);
dstate_x = jacobian(dstate, x_);
```

### Step 5: Export Functions
Save the derived models as MATLAB files. Keep filter models (`f`, `h`, `dstate`) separate from the truth models (`dstate_true`, `h_true`) for generality:
```matlab
matlabFunction(f,        'file', 'f',           'vars', {x_, u_, t, param});
matlabFunction(f_x,      'file', 'f_x',         'vars', {x_, u_, t, param});
matlabFunction(f_u,      'file', 'f_u',         'vars', {x_, u_, t, param});
matlabFunction(h,        'file', 'h',           'vars', {x_, u_, t, param});
matlabFunction(h_x,      'file', 'h_x',         'vars', {x_, u_, t, param});
matlabFunction(h_u,      'file', 'h_u',         'vars', {x_, u_, t, param});
matlabFunction(dstate,   'file', 'dstate',      'vars', {x_, u_, t, param});
matlabFunction(dstate_x, 'file', 'dstate_x',    'vars', {x_, u_, t, param});

% Separate true representations for truth simulation (can use different variables/parameters)
matlabFunction(h_true,      'file', 'h_true',      'vars', {x_true, u_true, t, param_true});
matlabFunction(dstate_true, 'file', 'dstate_true', 'vars', {x_true, u_true, t, param_true});
```

---

## 3. Writing `main_numeric.m`

The numeric script initializes parameters, configures noise standard deviations, executes the filter/smoother passes, and handles parameter estimation.

### Step 1: Initialize System Struct `S`
Define physical parameters, time steps, initial state estimates, and noise standard deviations:
```matlab
S.param = [1.0; 0.5; 9.81; 0.001]; % m, l, g, Delta_t
S.param_true = S.param;             % Parameter truth representation
S.Delta_t = 0.001;
S.t_end = 2.0;

% Initial conditions
S.t_0 = 0.0;
S.x_true_0 = [0.1; 0.0];            % True starting physical state
S.mu_x_0 = [0.15; 0.0];             % Filter starting state estimate
S.sigma_x_0 = [0.01; 0.01];         % Filter initial standard deviations
S.Sigma2_x_0 = diag(S.sigma_x_0.^2);

% Noise Standard Deviations (can be scalars or vectors)
S.sigma_w_x = [0.001; 0.001];       % Process noise std.
S.sigma_u = [];                     % Input noise std.
S.sigma_z = [0.01; 0.02];           % Sensor noise std.
S.sigma_v_x = [0.005; 0.005];       % Sensor process-coupling noise std.
```

### Step 2: Configure Logging & Run
Specify which fields should be saved to history in `datalogging_string`, then execute the forward pass and smoothers directly:
```matlab
datalogging_string = {'t'; 'x_true'; 'mu_x'; 'Sigma2_x'; 'u_meas'; 'z_meas'; 'sigma_x'};

% 1. Forward Filter Pass
Simulation = IEKF(S, datalogging_string);

% 2. Backward Smoother Pass (IEKS)
Simulation = IEKS(Simulation, S);

% 3. Iterated Backward Smoother Pass (IIEKS) - 4 iterations
% Simulation = IIEKS(Simulation, S, 4);
```

### Step 3: Unpack & Plot
Unpack the `Simulation` struct arrays to standard vectors/matrices for graphing:
```matlab
unpack_simulation(Simulation);

% Plotted variables are suffix-expanded to arrays (e.g. mu_x_series, x_true_series)
plot(t_series, x_true_series(:,1), 'k-', 'DisplayName', 'True');
hold on;
plot(t_series, mu_x_series(:,1), 'b--', 'DisplayName', 'Filter');
plot(t_series, mu_x_sm_series(:,1), 'r:', 'DisplayName', 'Smoother');
legend();
```

---

## 4. Parameter Estimation via Log-Likelihood

To optimize parameters (like noise standard deviations or model coefficients) using maximum likelihood estimation, define the parameter configuration inside the struct `S` and pass it to `logL_IEKF`.

### Configure Smoother Type in Likelihood
You can choose the state estimation trajectory used to compute the residuals inside `logL_IEKF`:
* **`S.smoother_in_likelihood = 'none'`**: Evaluates likelihood using the forward filter's predictions (standard prediction error).
* **`S.smoother_in_likelihood = 'IEKS'`**: Evaluates the joint energy function of the smoothed trajectory.
* **`S.smoother_in_likelihood = 'IIEKS'`**: Evaluates the joint energy function of the iterated smoothed trajectory.

### Run Optimization
```matlab
% Set up optimization variables (e.g. optimizing process noise sigma_w_x)
S.smoother_in_likelihood = 'none';

% Anonymous objective function that updates the struct and computes logL
obj_fun = @(theta) logL_IEKF(update_struct(S, 'sigma_w_x', theta));

% Optimize
initial_guess = [0.01; 0.01];
optimal_sigma = fminsearch(obj_fun, initial_guess);
```
