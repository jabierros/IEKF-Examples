function Simulation = IIEKS(Simulation, S, N_iter)
% IIEKS Performs trajectory-iterated backward pass and smoothing fusion.
%
% Simulation = IIEKS(Simulation, S, N_iter)

    if nargin < 3
        N_iter = 4;
    end

    N = Simulation.k;
    if N <= 0
        return;
    end
    
    n_x = size(S.mu_x_0, 1);
    
    % Reconstruct covariances from standard deviations or factors
    Sigma2_w_x = get_cov(S.sigma_w_x);
    Sigma2_u = get_cov(S.sigma_u);
    Sigma2_z = get_cov(S.sigma_z);
    Sigma2_v_x = get_cov(S.sigma_v_x);
    
    % Extract logged data cell arrays
    t = Simulation.data.t;
    u_meas = Simulation.data.u_meas;
    z_meas = Simulation.data.z_meas;
    
    % Initialize nominal smoothed trajectory with forward filter trajectory
    mu_x_sm = Simulation.data.mu_x;
    
    for iter = 1:N_iter
        % --- Nominal-Trajectory Forward Pass ---
        mu_x_fwd = cell(N, 1);
        Sigma2_x_fwd = cell(N, 1);
        mu_x_fwd{1} = S.mu_x_0;
        Sigma2_x_fwd{1} = S.Sigma2_x_0;
        
        for k = 1:N-1
            % Prediction from k to k+1
            u_meas_pred = u_meas{k+1};
            t_pred = t{k};
            
            f_x_ = f_x(mu_x_sm{k}, u_meas_pred, t_pred, S.param);
            f_u_ = f_u(mu_x_sm{k}, u_meas_pred, t_pred, S.param);
            f_val = f(mu_x_sm{k}, u_meas_pred, t_pred, S.param);
            
            Q = f_u_ * Sigma2_u * f_u_' + Sigma2_w_x;
            
            % Taylor linearized mean and covariance propagation
            mu_x_pred = f_val + f_x_ * (mu_x_fwd{k} - mu_x_sm{k});
            Sigma2_x_pred = f_x_ * Sigma2_x_fwd{k} * f_x_' + Q;
            
            I_x_pred = inv(Sigma2_x_pred);
            i_x_pred = I_x_pred * mu_x_pred;
            
            % Update at k+1
            z_meas_upd = z_meas{k+1};
            t_upd = t{k+1};
            
            h_x_ = h_x(mu_x_sm{k+1}, u_meas_pred, t_upd, S.param);
            h_u_ = h_u(mu_x_sm{k+1}, u_meas_pred, t_upd, S.param);
            h_val = h(mu_x_sm{k+1}, u_meas_pred, t_upd, S.param);
            
            R = Sigma2_z + h_u_ * Sigma2_u * h_u_' + Sigma2_v_x;
            inv_R = inv(R);
            
            I_x_obs = h_x_' * inv_R * h_x_;
            i_x_obs = h_x_' * inv_R * (z_meas_upd - (h_val - h_x_ * mu_x_sm{k+1}));
            
            I_x = I_x_pred + I_x_obs;
            i_x = i_x_pred + i_x_obs;
            
            mu_x_fwd{k+1} = inv(I_x) * i_x;
            Sigma2_x_fwd{k+1} = inv(I_x);
        end
        
        % --- Nominal-Trajectory Backward Pass & Fusion ---
        I_x_pred_bck = zeros(n_x, n_x);
        i_x_pred_bck = zeros(n_x, 1);
        
        for k = N:-1:1
            mu_x_nominal_k = mu_x_sm{k};
            if k > 1
                mu_x_nominal_prev = mu_x_sm{k-1};
                u_meas_prev = u_meas{k};
                t_prev = t{k-1};
            else
                mu_x_nominal_prev = [];
                u_meas_prev = [];
                t_prev = [];
            end
            
            [mu_x_sm{k}, Sigma2_x_sm{k}, I_x_pred_bck, i_x_pred_bck] = IEKS_step(...
                k, mu_x_nominal_k, mu_x_nominal_prev, mu_x_fwd{k}, Sigma2_x_fwd{k}, ...
                u_meas{k}, z_meas{k}, t{k}, u_meas_prev, t_prev, I_x_pred_bck, i_x_pred_bck, S);
        end
    end
    
    % Store final iterated smoothed results in Simulation data
    Simulation.data.mu_x_sm = mu_x_sm;
    Simulation.data.Sigma2_x_sm = Sigma2_x_sm;
    
    % Append names for unpack_simulation to dynamically pick up
    if ~any(strcmp(Simulation.names, 'mu_x_sm'))
        Simulation.names{end+1} = 'mu_x_sm';
    end
    if ~any(strcmp(Simulation.names, 'Sigma2_x_sm'))
        Simulation.names{end+1} = 'Sigma2_x_sm';
    end
end

function cov_matrix = get_cov(sigma)
    if size(sigma, 2) == 1
        cov_matrix = diag(sigma.^2);
    else
        cov_matrix = sigma * sigma';
    end
end
