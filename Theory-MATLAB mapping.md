# Theory-MATLAB Mapping Document

This document defines the 1-to-1 correspondence between the mathematical concepts in the Kalman Filter / Information Extended Kalman Filter (IEKF) / Information Extended Kalman Smoother (IEKS) theory and the variables used in the MATLAB codebase.

## Variable Mapping Table

| **Concept** | **Theory / Math** | **MATLAB Variable (Current Step)** | **MATLAB Array (History)** |
| :--- | :--- | :--- | :--- |
| **Time & Params** | $t_k$, $\Delta t$, $\mathbf{p}$ | `t`, `Delta_t`, `param` | `t_series` |
| **True Reality** | $\mathbf{x}^{tr}$, $\mathbf{z}^{tr}$, $\mathbf{u}^{tr}$ | `x_true`, `z_true`, `u_true` | `x_true_series`, `u_true_series`, `z_true_series` |
| **Observations** | $\mathbf{z}$, $\mathbf{u}$ | `z_meas`, `u_meas` | `z_meas_series`, `u_meas_series` |
| **Jacobians (Evaluated)** | $\mathbf{f_x}$, $\mathbf{f_u}$, $\mathbf{h_x}$, $\mathbf{h_u}$ | `f_x_`, `f_u_`, `h_x_`, `h_u_` | N/A |
| **Noise Std. Devs** | $\sigma_{\mathbf{w}}, \sigma_{\mathbf{v}}, \sigma_{\mathbf{u}}, \sigma_{\mathbf{z}}$ | `sigma_w_x`, `sigma_v_x`, `sigma_u`, `sigma_z` | N/A |
| **Projected Input Noise** | $\boldsymbol{\Sigma}_{\mathbf{f_u u}}^2$, $\boldsymbol{\Sigma}_{\mathbf{h_u u}}^2$ | `Sigma2_f_u_u`, `Sigma2_h_u_u` | N/A |
| **Prediction (Source 1)** | $\hat{\boldsymbol{\mu}}_{\mathbf{x}}^{pred}$, $\mathbf{I}_{\mathbf{x}}^{pred}$, $\hat{\mathbf{i}}_{\mathbf{x}}^{pred}$ | `mu_x_pred`, `I_x_pred`, `i_x_pred` | N/A |
| **Observation (Source 2)** | $\mathbf{I}_{\mathbf{x}}^{obs}$, $\hat{\mathbf{i}}_{\mathbf{x}}^{obs}$ | `I_x_obs`, `i_x_obs` | N/A |
| **Aggregated Updated** | $\hat{\boldsymbol{\mu}}_{\mathbf{x}}$, $\boldsymbol{\Sigma}_{\mathbf{x}}^2$, $\mathbf{I}_{\mathbf{x}}$, $\hat{\mathbf{i}}_{\mathbf{x}}$ | `mu_x`, `Sigma2_x`, `I_x`, `i_x` | `mu_x_series`, `Sigma2_x_series`, `I_x_series` |
| **Backward Prediction** | $\mathbf{I}_{\mathbf{x}}^{pred,bck}$, $\hat{\mathbf{i}}_{\mathbf{x}}^{pred,bck}$ | `I_x_pred_bck`, `i_x_pred_bck` | N/A |
| **Backward Updated** | $\mathbf{I}_{\mathbf{x}}^{bck}$, $\hat{\mathbf{i}}_{\mathbf{x}}^{bck}$ | `I_x_bck`, `i_x_bck` | N/A |
| **Smoothed Fusion** | $\hat{\boldsymbol{\mu}}_{\mathbf{x}}^{sm}$, $\boldsymbol{\Sigma}_{\mathbf{x}}^{2,sm}$, $\mathbf{I}_{\mathbf{x}}^{sm}$, $\hat{\mathbf{i}}_{\mathbf{x}}^{sm}$ | `mu_x_sm`, `Sigma2_x_sm`, `I_x_sm`, `i_x_sm` | `mu_x_sm_series`, `Sigma2_x_sm_series`, `I_x_sm_series` |
| **Statistical Analysis** | $RMS(\epsilon_{\mathbf{x}})$ | `sqrt_lim_mu_x_error_squared_mean` | `mu_x_error_series` |

## Notation Details & Conceptual Separations

### 1. Symbolic Notation (`x_`, `u_`, `z_`)
* In the symbolic script [main_symbolic_EKF.m](file:///home/jros/Sync/Kalman/IEKF-Examples/Double_Pendulum/main_symbolic_EKF.m), `x_` is used for the symbolic state vector (instead of `x`, to avoid naming conflicts with spatial coordinates).
* To be fully coherent, the symbolic input vector is named `u_` and exported functions are generated using the variables list `{x_, u_, t, param}` or `{x_true, u_true, t, param_true}`.

### 2. True vs. Filter Continuous Models
We conceptually split continuous system modeling:
* `dstate` represents the continuous system derivative model used in the filter (e.g. for discretization in the prediction step). It depends on `param`.
* `dstate_true` represents the true system derivative model used to simulate reality (integrated via `ode45`). It depends on `param_true` and uses `x_true` and `u_true` symbolic coordinates (which are coincident with `x_` and `u_` in the symbolic derivation).

### 3. Measurements Prefix (`_meas`)
* Variables representing observations/measurements are suffixed with `_meas` (e.g. `u_meas`, `z_meas`, `u_meas_series`, `z_meas_series`) to clearly delineate them from states or true physical properties.

### 4. Two-Filter Smoothing Notation (`_bck` and `_sm`)
* Variables associated with the backward filter pass are suffixed with `_bck` (for updated) or `_pred_bck` (for predicted back in time).
* Smoothed variables resulting from the fusion of the forward updated and backward predicted information are suffixed with `_sm` (e.g., `mu_x_sm`, `Sigma2_x_sm`).

### 5. Iterated Smoothing (IIEKS)
* Unlike the standard EKS/IEKS which linearizes once along the forward trajectory, the Iterated Information Extended Kalman Smoother (IIEKS) performs multiple forward-backward optimization passes.
* In each iteration, both the nominal forward pass and the backward pass equations are linearized (evaluating process and measurement Jacobians) around the *latest smoothed trajectory* $\hat{\boldsymbol{\mu}}_{\mathbf{x}}^{sm}$ from the previous iteration. This centers the Taylor series expansions closer to the true states, reducing linearization error and improving the marginal log-likelihood.

### 6. Architectural Structure & Likelihood Option
* **Pass-level vs. Step-level Functions:**
  * `IEKF.m` and `IEKS.m` are **pass-level** functions that run the entire trajectory forward or backward and return/store history in the `Simulation` struct.
  * `IEKF_step.m` and `IEKS_step.m` are **step-level** functions that run the filtering/smoothing equations for a single time step $k$. They encapsulate the mathematical formulas (prediction, update, fusion) cleanly.
* **Log-Likelihood Smoother Options:**
  * The negative log-likelihood/joint energy can be computed with various settings using `S.smoother_in_likelihood`:
    * `'none'`: Forward prediction error loop using `IEKF_step(S)` (standard log-likelihood).
    * `'IEKS'`: Full pass via `IEKF`, smoothed pass via `IEKS`, and evaluation of the joint state-measurement negative log-likelihood (energy) of the smoothed trajectory.
    * `'IIEKS'`: Full pass via `IEKF`, iterated smoothed pass via `IIEKS`, and evaluation of the joint negative log-likelihood (energy) of the iterated smoothed trajectory.

