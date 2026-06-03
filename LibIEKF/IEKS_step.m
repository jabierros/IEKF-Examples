function [mu_x_sm_k, Sigma2_x_sm_k, I_x_pred_bck_prev, i_x_pred_bck_prev] = IEKS_step(...
    k, mu_x_nominal_k, mu_x_nominal_prev, mu_x_fwd_k, Sigma2_x_fwd_k, ...
    u_meas_k, z_meas_k, t_k, u_meas_prev, t_prev, I_x_pred_bck, i_x_pred_bck, S)
% IEKS_STEP Performs one step of the backward pass and smoothing fusion.
%
% [mu_x_sm_k, Sigma2_x_sm_k, I_x_pred_bck_prev, i_x_pred_bck_prev] = IEKS_step(...)

    % Reconstruct covariances
    Sigma2_w_x = get_cov(S.sigma_w_x);
    Sigma2_u = get_cov(S.sigma_u);
    Sigma2_z = get_cov(S.sigma_z);
    Sigma2_v_x = get_cov(S.sigma_v_x);
    
    n_x = size(mu_x_nominal_k, 1);
    
    % 1. Backward Update at step k
    h_x_val = h_x(mu_x_nominal_k, u_meas_k, t_k, S.param);
    h_u_val = h_u(mu_x_nominal_k, u_meas_k, t_k, S.param);
    h_ = h(mu_x_nominal_k, u_meas_k, t_k, S.param);
    
    R = Sigma2_z + h_u_val * Sigma2_u * h_u_val' + Sigma2_v_x;
    inv_R = inv(R);
    
    I_x_obs = h_x_val' * inv_R * h_x_val;
    i_x_obs = h_x_val' * inv_R * (z_meas_k - (h_ - h_x_val * mu_x_nominal_k));
    
    I_x_bck = I_x_pred_bck + I_x_obs;
    i_x_bck = i_x_pred_bck + i_x_obs;
    
    % 2. Smoothing Fusion at step k
    I_x_fwd = inv(Sigma2_x_fwd_k);
    i_x_fwd = I_x_fwd * mu_x_fwd_k;
    
    I_x_sm = I_x_fwd + I_x_pred_bck;
    i_x_sm = i_x_fwd + i_x_pred_bck;
    
    mu_x_sm_k = inv(I_x_sm) * i_x_sm;
    Sigma2_x_sm_k = inv(I_x_sm);
    
    % 3. Backward Prediction from k back to k-1
    if k > 1
        f_x_val = f_x(mu_x_nominal_prev, u_meas_prev, t_prev, S.param);
        f_u_val = f_u(mu_x_nominal_prev, u_meas_prev, t_prev, S.param);
        f_val = f(mu_x_nominal_prev, u_meas_prev, t_prev, S.param);
        
        Q = f_u_val * Sigma2_u * f_u_val' + Sigma2_w_x;
        inv_Q = inv(Q);
        
        M = inv(I_x_bck + inv_Q);
        W = I_x_bck - I_x_bck * M * I_x_bck;
        
        I_x_pred_bck_prev = f_x_val' * W * f_x_val;
        f_res = f_val - f_x_val * mu_x_nominal_prev;
        i_x_pred_bck_prev = f_x_val' * ((eye(n_x) - I_x_bck * M) * i_x_bck - W * f_res);
    else
        I_x_pred_bck_prev = zeros(n_x, n_x);
        i_x_pred_bck_prev = zeros(n_x, 1);
    end
end

function cov_matrix = get_cov(sigma)
    if size(sigma, 2) == 1
        cov_matrix = diag(sigma.^2);
    else
        cov_matrix = sigma * sigma';
    end
end
