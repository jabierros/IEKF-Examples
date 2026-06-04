function [mu_x_next, Sigma2_x_next, mu_x_pred, Sigma2_x_pred] = IEKF_step(...
    mu_x, Sigma2_x, u_meas_prev, u_meas, z_meas, t_prev, t, Delta_t, ...
    Sigma_w, Sigma_v, Sigma_u, Sigma_z, param)
% IEKF_STEP Performs a single prediction and update step of the Iterated Extended Kalman Filter.
%
% [mu_x_next, Sigma2_x_next, mu_x_pred, Sigma2_x_pred] = IEKF_step(...)

    % Reconstruct covariances from Cholesky factors
    Sigma2_w = Sigma_w * Sigma_w';
    Sigma2_u = Sigma_u * Sigma_u';
    Sigma2_z = Sigma_z * Sigma_z';
    Sigma2_v = Sigma_v * Sigma_v';

    % 1. Prediction (from t_prev to t, using control input from t_prev)
    f_u = f_u_(mu_x, u_meas_prev, t_prev, param, Delta_t);
    f_x = f_x_(mu_x, u_meas_prev, t_prev, param, Delta_t);
    f = f_(mu_x, u_meas_prev, t_prev, param, Delta_t);

    Q = f_u * Sigma2_u * f_u' + Sigma2_w;

    mu_x_pred = f;
    Sigma2_x_pred = f_x * Sigma2_x * f_x' + Q;
    I_x_pred = inv(Sigma2_x_pred);
    i_x_pred = I_x_pred * mu_x_pred;

    % 2. Projection / Observation (at t, using control input from t)
    h_u = h_u_(mu_x_pred, u_meas, t, param);
    h_x = h_x_(mu_x_pred, u_meas, t, param);
    h = h_(mu_x_pred, u_meas, t, param);

    R = Sigma2_z + h_u * Sigma2_u * h_u' + Sigma2_v;
    inv_R = inv(R);

    I_x_obs = h_x' * inv_R * h_x;
    i_x_obs = h_x' * inv_R * (z_meas - (h - h_x * mu_x_pred));

    % 3. Fusion / Information Weighted Average
    I_x_next = I_x_pred + I_x_obs;
    i_x_next = i_x_pred + i_x_obs;

    Sigma2_x_next = inv(I_x_next);
    mu_x_next = Sigma2_x_next * i_x_next;
end
