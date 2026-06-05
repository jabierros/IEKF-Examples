# LibIEKF User Manual: Writing Symbolic & Numeric Scripts

This manual explains how to write new filter/smoother applications using the **Information Extended Kalman Filter and Smoother (`LibIEKF`)** library. 

---

## 1. Library Architecture & Requirements

`LibIEKF` uses an in-memory, struct-based architecture where parameters, state representations, and simulation options are separated into three structures: `KF`, `TrueSystem`, and `SimOpts`.

### Required Template Functions
Every system directory (e.g., `Double_Pendulum`) uses the following interface functions located under `LibIEKF/Template` to supply inputs and measurements:
1. **`get_u.m`**: Fetches or generates the command input `u_meas`.
   ```matlab
   function [u_meas, TrueSystem, SimOpts] = get_u(TrueSystem, SimOpts)
   ```
2. **`get_z.m`**: Fetches or generates the measurement vector `z_meas`.
   ```matlab
   function [z_meas, TrueSystem, SimOpts] = get_z(TrueSystem, SimOpts)
   ```
3. **`get_x_true.m`**: Integrates continuous physics equations using the true process model (`dstate_true_.m`) to generate the true physical trajectory `x_true`.
   ```matlab
   function [x_true, TrueSystem, SimOpts] = get_x_true(TrueSystem, SimOpts)
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
param = [m; l; g];          % System parameter vector (excludes Delta_t)
param_true = param;
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
Save derived models with trailing underscores (`_`) to prevent name shadowing. Keep filter models separate from truth models:
```matlab
matlabFunction(f,        'file', 'f_',           'vars', {x_, u_, t, param, Delta_t});
matlabFunction(f_x,      'file', 'f_x_',         'vars', {x_, u_, t, param, Delta_t});
matlabFunction(f_u,      'file', 'f_u_',         'vars', {x_, u_, t, param, Delta_t});
matlabFunction(h,        'file', 'h_',           'vars', {x_, u_, t, param});
matlabFunction(h_x,      'file', 'h_x_',         'vars', {x_, u_, t, param});
matlabFunction(h_u,      'file', 'h_u_',         'vars', {x_, u_, t, param});
matlabFunction(dstate,   'file', 'dstate_',      'vars', {x_, u_, t, param});
matlabFunction(dstate_x, 'file', 'dstate_x_',    'vars', {x_, u_, t, param});

% Separate true representations for truth simulation (with symbolic variables u_true and x_true)
x_true = x_;
syms u_true real
dstate_true = subs(dstate, u_, u_true);
h_true = subs(h, u_, u_true);

matlabFunction(h_true,      'file', 'h_true_',      'vars', {x_true, u_true, t, param_true});
matlabFunction(dstate_true, 'file', 'dstate_true_', 'vars', {x_true, u_true, t, param_true});
matlabFunction(sym(0),      'file', 'u_true_func_', 'vars', {t});
```

---

## 3. Writing `main_numeric.m`

The numeric script initializes parameters, configures noise standard deviations, executes the filter/smoother passes, and handles parameter estimation.

### Step 1: Initialize Structs
Define physical parameters, time steps, initial state estimates, and noise standard deviations in three decoupled structs:
```matlab
% 1. Simulation Options
SimOpts = struct();
SimOpts.t_0 = 0.0;
SimOpts.t_end = 2.0;
SimOpts.t = 0.0;
SimOpts.t_prev = 0.0;

% 2. Filter Configuration
KF = struct();
KF.Delta_t = 0.001;
KF.param = [1.0; 0.5; 9.81];        % m, l, g
KF.mu_x_0 = [0.15; 0.0];             % Filter starting state estimate
KF.sigma_x_0 = [0.01; 0.01];         % Filter initial standard deviations
KF.Sigma_w = diag([0.001; 0.001]);   % Process noise Cholesky matrix
KF.Sigma_v = diag([0.005; 0.005]);   % Sensor process-coupling Cholesky matrix
KF.Sigma_u = diag([]);               % Input noise Cholesky matrix
KF.Sigma_z = diag([0.01; 0.02]);     % Sensor noise Cholesky matrix

% 3. True System configuration
TrueSystem = struct();
TrueSystem.param_true = KF.param;    % True parameters
TrueSystem.x_true_0 = [0.1; 0.0];    % True starting state
TrueSystem.x_true = [0.1; 0.0];
TrueSystem.x_true_prev = [0.1; 0.0];
TrueSystem.sigma_u_true = [];        % True input noise level
TrueSystem.sigma_z_true = [0.01; 0.02]; % True sensor noise level
```

### Step 2: Configure Logging & Run
Specify which fields should be saved to history in `datalogging_string`, then execute the forward pass and smoothers:
```matlab
datalogging_string = {'t'; 'x_true'; 'mu_x'; 'Sigma2_x'; 'u_meas'; 'z_meas'; 'sigma_x'};

% 1. Forward Filter Pass
FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);

% 2. Backward Smoother Pass (IEKS)
FilterResults = IEKS(FilterResults, KF, SimOpts);

% 3. Iterated Backward Smoother Pass (IIEKS) - 4 iterations
% FilterResults = IIEKS(FilterResults, KF, SimOpts, 4);
```

### Step 3: Unpack & Plot
Unpack the `FilterResults` struct arrays to standard vectors/matrices for graphing:
```matlab
unpack_simulation(FilterResults);

% Plotted variables are suffix-expanded to arrays (e.g. mu_x_series, x_true_series)
plot(t_series, x_true_series(:,1), 'k-', 'DisplayName', 'True');
hold on;
plot(t_series, mu_x_series(:,1), 'b--', 'DisplayName', 'Filter');
plot(t_series, mu_x_sm_series(:,1), 'r:', 'DisplayName', 'Smoother');
legend();
```

---

## 4. Parameter Estimation via Log-Likelihood

To optimize parameters (like noise standard deviations or model coefficients) using maximum likelihood estimation, pass updated structures to `logL_IEKF`:

```matlab
% Set up optimization variables (e.g. optimizing process noise sigma_w)
theta = [sigma_w];

% Anonymous objective function that updates the struct and computes logL
obj_fun = @(theta) logL_IEKF(update_struct(KF, 'sigma_w', theta), TrueSystem, SimOpts);

% Optimize
options = optimset('PlotFcns', @optimplotfval);
optimal_sigma = fminsearch(obj_fun, theta, options);
```
