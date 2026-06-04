function [mu_x_sm_k, Sigma2_x_sm_k, I_x_pred_bck_prev, i_x_pred_bck_prev] = IEKS_step(...
    k, mu_x_nominal_k, mu_x_nominal_prev, mu_x_fwd_k, Sigma2_x_fwd_k, ...
    u_meas_k, z_meas_k, t_k, u_meas_prev, t_prev, I_x_pred_bck, i_x_pred_bck, ...
    Delta_t, Sigma_w, Sigma_v, Sigma_u, Sigma_z, param)
% IEKS_STEP Performs one step of the backward pass and smoothing fusion.
%
% [mu_x_sm_k, Sigma2_x_sm_k, I_x_pred_bck_prev, i_x_pred_bck_prev] = IEKS_step(...)

    % Reconstruct covariances
    Sigma2_w = Sigma_w * Sigma_w';
    Sigma2_u = Sigma_u * Sigma_u';
    Sigma2_z = Sigma_z * Sigma_z';
    Sigma2_v = Sigma_v * Sigma_v';
    
    n_x = size(mu_x_nominal_k, 1);
    
    % 1. Backward Update at step k
    h_x = h_x_(mu_x_nominal_k, u_meas_k, t_k, param);
    h_u = h_u_(mu_x_nominal_k, u_meas_k, t_k, param);
    h = h_(mu_x_nominal_k, u_meas_k, t_k, param);
    
    R = Sigma2_z + h_u * Sigma2_u * h_u' + Sigma2_v;
    inv_R = inv(R);
    
    I_x_obs = h_x' * inv_R * h_x;
    i_x_obs = h_x' * inv_R * (z_meas_k - (h - h_x * mu_x_nominal_k));
    
    I_x_bck = I_x_pred_bck + I_x_obs;
    i_x_bck = i_x_pred_bck + i_x_obs;
    
    % 2. Smoothing Fusion at step k
    I_x_fwd = inv(Sigma2_x_fwd_k);
    i_x_fwd = I_x_fwd * mu_x_fwd_k;
    
    I_x_sm = I_x_fwd + I_x_pred_bck;
    i_x_sm = i_x_fwd + i_x_pred_bck;
    
    Sigma2_x_sm_k = inv(I_x_sm);
    mu_x_sm_k = Sigma2_x_sm_k * i_x_sm;
    
    % 3. Backward Prediction from k back to k-1
    if k > 1
        f_x = f_x_(mu_x_nominal_prev, u_meas_prev, t_prev, param, Delta_t);
        f_u = f_u_(mu_x_nominal_prev, u_meas_prev, t_prev, param, Delta_t);
        f = f_(mu_x_nominal_prev, u_meas_prev, t_prev, param, Delta_t);
        
        Q = f_u * Sigma2_u * f_u' + Sigma2_w;
        inv_Q = inv(Q);
        
        M = inv(I_x_bck + inv_Q);
        W = I_x_bck - I_x_bck * M * I_x_bck;
        
        I_x_pred_bck_prev = f_x' * W * f_x;
        f_res = f - f_x * mu_x_nominal_prev;
        i_x_pred_bck_prev = f_x' * ((eye(n_x) - I_x_bck * M) * i_x_bck - W * f_res);
    else
        I_x_pred_bck_prev = zeros(n_x, n_x);
        i_x_pred_bck_prev = zeros(n_x, 1);
    end
end
