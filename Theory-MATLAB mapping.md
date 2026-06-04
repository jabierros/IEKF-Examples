# Theory-MATLAB Notation Guideline

This document defines the mathematical and programming conventions for variables, files, and functions in the Kalman Filter / Information Extended Kalman Filter (IEKF) / Information Extended Kalman Smoother (IEKS) workspace. Rather than a rigid dictionary, it serves as a compositional style guide to maintain clarity, avoid variable conflicts, and preserve exact mapping to the filter theory equations.

---

## 1. Core Mathematical Root Names
The primary variables in theory map to lowercase roots in MATLAB:
* State: $\mathbf{x} \to$ `x`
* Measurement / Output: $\mathbf{z} \to$ `z`
* Control / Input: $\mathbf{u} \to$ `u`
* Process Transition Function: $\mathbf{f} \to$ `f`
* Measurement Function: $\mathbf{h} \to$ `h`
* Continuous Derivative Function: $\dot{\mathbf{x}} \to$ `dstate`

---

## 2. Capitalization for Dimensionality (Matrices vs. Vectors/Scalars)
Casing distinguishes vectors or scalars from matrices and operators:
* **Lowercase (`sigma`, `mu`, `i`):** Used for vectors or scalars.
  * E.g., `sigma_w` is the standard deviation *vector* of the process noise.
  * E.g., `mu_x` is the state mean *vector* $\boldsymbol{\mu}_{\mathbf{x}}$.
  * E.g., `i_x` is the information *vector* $\mathbf{i}_{\mathbf{x}}$.
* **Capitalized (`Sigma`, `I`):** Used for matrices.
  * E.g., `Sigma_w` is the Cholesky *matrix* factor $\boldsymbol{\Sigma}_{\mathbf{w}}$.
  * E.g., `Sigma2_w` is the covariance *matrix* $\boldsymbol{\Sigma}_{\mathbf{w}}^2 = \boldsymbol{\Sigma}_{\mathbf{w}}\boldsymbol{\Sigma}_{\mathbf{w}}^\top$.
  * E.g., `I_x` is the Information *matrix* $\mathbf{I}_{\mathbf{x}}$.

---

## 3. Compositional Suffix/Modifier System
Variable names are constructed by appending suffixes to core roots:
$$\text{MATLAB Name} = \langle\text{Root}\rangle\_[\text{State/Type}]\_[\text{Step/Epoch}]\_[\text{History}]$$

### State / Type Suffixes
* `_true`: True physical value (in simulation).
  * *Note on Epistemology:* In the physical simulation, the actual state/measurement/input representing reality are named `x_true`, `z_true`, `u_true`. These represent the true reality and are inaccessible to the filter. Functions and parameters without `_true` represent the internal model of the filter, allowing simulation of reality-model discrepancies (e.g. model mismatch).
* `_meas`: Sensor readings / measurements (accessible to the filter), e.g. `z_meas`, `u_meas`.
* `_pred`: Predicted value (prior), e.g., `mu_x_pred`, `Sigma2_x_pred`.
* `_obs`: Measurement update component (observation info), e.g., `I_x_obs`, `i_x_obs`.
* `_bck`: Backward pass, e.g., `I_x_bck`, `i_x_bck`.
* `_pred_bck`: Backward prediction, e.g., `I_x_pred_bck`, `i_x_pred_bck`.
* `_sm`: Smoothed value, e.g., `mu_x_sm`, `Sigma2_x_sm`.
* `_error`: Statistical discrepancy or error, e.g., `mu_x_error`.

### Step / Epoch Suffixes
Sequence indexes in time loops are denoted with subindexes:
* `_k`: Step $k$, e.g., `u_meas_k`, `mu_x_k`.
* `_kp1`: Step $k+1$, e.g., `u_meas_kp1`, `mu_x_kp1`.
* `_km1`: Step $k-1$, e.g., `u_meas_km1`, `mu_x_km1`.
* `_prev` and `_next`: Local variable overrides used strictly inside step-level functions (like `IEKF_step.m` or `IEKS_step.m`) where there is no loop index `k` and the inputs/outputs must be disambiguated.

### History Suffixes
* `_series`: A logged time-series collection of a variable, e.g., `mu_x_series`, `mu_x_error_series`.

---

## 4. Function File Names vs. Evaluated Variables
To avoid MATLAB naming conflicts when calling a function and storing its return value, we append a trailing underscore (`_`) to the **exported function file names**.
This allows the local variables holding their evaluated values to match the exact mathematical symbols:

| Mathematical Concept | Function File | Evaluation in Code |
| :--- | :--- | :--- |
| Process Transition $\mathbf{f}$ | `f_.m` | `f = f_(x, u, t, param, Delta_t)` |
| Process Jacobian $\mathbf{f_x}$ | `f_x_.m` | `f_x = f_x_(x, u, t, param, Delta_t)` |
| Process Control Jacobian $\mathbf{f_u}$ | `f_u_.m` | `f_u = f_u_(x, u, t, param, Delta_t)` |
| Measurement Function $\mathbf{h}$ | `h_.m` | `h = h_(x, u, t, param)` |
| Measurement Jacobian $\mathbf{h_x}$ | `h_x_.m` | `h_x = h_x_(x, u, t, param)` |
| Measurement Control Jacobian $\mathbf{h_u}$ | `h_u_.m` | `h_u = h_u_(x, u, t, param)` |
| Continuous Derivative $\dot{\mathbf{x}}$ | `dstate_` | `dstate = dstate_(x, u, t, param)` |
| Continuous Jacobian $\dot{\mathbf{x}}_{\mathbf{x}}$ | `dstate_x_` | `dstate_x = dstate_x_(x, u, t, param)` |
| True Reality Derivative | `dstate_true_` | `dstate_true = dstate_true_(x, u, t, param)` |
| True Reality Measurement | `h_true_` | `h_true = h_true_(x, u, t, param)` |
| True Input Function | `u_true_func_` | `u_true = u_true_func_(t)` |

*Note:* Evaluated variables can still be combined with step-based modifiers as needed, e.g., `f_x_k = f_x_(mu_x_k, ...)` or `f_pred = f_(...)`.

---

## 5. Architectural Conventions
* **Pass-Level vs. Step-Level Separation:**
  * `IEKF.m` and `IEKS.m` are **pass managers** that handle state arrays over a trajectory and return logged results.
  * `IEKF_step.m` and `IEKS_step.m` are **step functions** evaluating the filter/smoother math on flat vectors/matrices at a single instant in time.
* **Decoupled Structs:**
  * `KF`: Filter parameters (`param`, `Delta_t`, noise Cholesky factors like `Sigma_w`).
  * `TrueSystem`: True physical reality settings (`param_true`, `x_true_0`, noise std `sigma_u_true`, `sigma_z_true`).
  * `SimOpts`: Time tracker (`t_0`, `t_end`, `t`, `t_prev`).
