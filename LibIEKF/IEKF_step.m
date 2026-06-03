function S = IEKF_step(S)
% IEKF_STEP Performs a single prediction and update step of the Iterated Extended Kalman Filter.
%
% S = IEKF_step(S)

    % Prediction
    [u_meas, S] = get_u(S);
    S.u_meas = u_meas;
    
    % Reconstruct covariances from standard deviations or factors
    Sigma2_w_x = get_cov(S.sigma_w_x);
    Sigma2_u = get_cov(S.sigma_u);
    Sigma2_z = get_cov(S.sigma_z);
    Sigma2_v_x = get_cov(S.sigma_v_x);
    
    % Process noise covariance (input noise + process noise sources)
    f_u_ = f_u(S.mu_x, S.u_meas, S.t, S.param);
    S.Q = f_u_ * Sigma2_u * f_u_' + Sigma2_w_x;
    
    f_x_ = f_x(S.mu_x, S.u_meas, S.t, S.param);
    
    S.mu_x_pred = f(S.mu_x, S.u_meas, S.t, S.param);
    S.Sigma2_x_pred = f_x_ * S.Sigma2_x * f_x_' + S.Q;
    I_x_pred = inv(S.Sigma2_x_pred);
    i_x_pred = I_x_pred * S.mu_x_pred;
    
    % Read Input and Sensor (k+1)
    S.t = S.t + S.Delta_t;
    
    [z_meas, S] = get_z(S);
    S.z_meas = z_meas;
    
    % Projection / Observation
    h_u_ = h_u(S.mu_x_pred, S.u_meas, S.t, S.param);
    
    % Measurement noise covariance (sensor noise + input noise + other noise sources)
    S.R = Sigma2_z + h_u_ * Sigma2_u * h_u_' + Sigma2_v_x;
    
    h_x_ = h_x(S.mu_x_pred, S.u_meas, S.t, S.param);
    i_x_obs = h_x_' * inv(S.R) * (S.z_meas - (h(S.mu_x_pred, S.u_meas, S.t, S.param) - h_x_ * S.mu_x_pred));
    I_x_obs = h_x_' * inv(S.R) * h_x_;
    
    % Fusion / Information Weighted Average
    i_x = i_x_pred + i_x_obs;
    I_x = I_x_pred + I_x_obs;
    S.mu_x = inv(I_x) * i_x;
    S.Sigma2_x = inv(I_x);
end

function cov_matrix = get_cov(sigma)
    if size(sigma, 2) == 1
        cov_matrix = diag(sigma.^2);
    else
        cov_matrix = sigma * sigma';
    end
end
