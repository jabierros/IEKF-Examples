function FilterResults = IIEKS(FilterResults, KF, SimOpts, N_iter)
% IIEKS Performs trajectory-iterated backward pass and smoothing fusion.
%
% FilterResults = IIEKS(FilterResults, KF, SimOpts, N_iter)

    if nargin < 4
        N_iter = 4;
    end

    N = FilterResults.k;
    if N <= 0
        return;
    end
    
    n_x = size(KF.mu_x_0, 1);
    
    % Reconstruct covariances from Cholesky factors
    Sigma2_w = KF.Sigma_w * KF.Sigma_w';
    Sigma2_u = KF.Sigma_u * KF.Sigma_u';
    Sigma2_z = KF.Sigma_z * KF.Sigma_z';
    Sigma2_v = KF.Sigma_v * KF.Sigma_v';
    
    % Extract logged data cell arrays
    t = FilterResults.data.t;
    u_meas = FilterResults.data.u_meas;
    z_meas = FilterResults.data.z_meas;
    
    % Initialize nominal smoothed trajectory with forward filter trajectory
    mu_x_sm = FilterResults.data.mu_x;
    
    % Set initial covariance
    if isvector(KF.sigma_x_0)
        Sigma2_x_0 = diag(KF.sigma_x_0.^2);
    else
        Sigma2_x_0 = KF.sigma_x_0 * KF.sigma_x_0';
    end
    
    for iter = 1:N_iter
        % --- Nominal-Trajectory Forward Pass ---
        mu_x_fwd = cell(N, 1);
        Sigma2_x_fwd = cell(N, 1);
        mu_x_fwd{1} = KF.mu_x_0;
        Sigma2_x_fwd{1} = Sigma2_x_0;
        
        for k = 1:N-1
            % Prediction from k to k+1 (using input from time t_k)
            u_meas_pred = u_meas{k};
            t_pred = t{k};
            
            f_x = f_x_(mu_x_sm{k}, u_meas_pred, t_pred, KF.param, KF.Delta_t);
            f_u = f_u_(mu_x_sm{k}, u_meas_pred, t_pred, KF.param, KF.Delta_t);
            f = f_(mu_x_sm{k}, u_meas_pred, t_pred, KF.param, KF.Delta_t);
            
            Q = f_u * Sigma2_u * f_u' + Sigma2_w;
            
            % Taylor linearized mean and covariance propagation
            mu_x_pred = f + f_x * (mu_x_fwd{k} - mu_x_sm{k});
            Sigma2_x_pred = f_x * Sigma2_x_fwd{k} * f_x' + Q;
            
            I_x_pred = inv(Sigma2_x_pred);
            i_x_pred = I_x_pred * mu_x_pred;
            
            % Update at k+1 (using input and measurement from time t_k+1)
            z_meas_upd = z_meas{k+1};
            u_meas_upd = u_meas{k+1};
            t_upd = t{k+1};
            
            h_x = h_x_(mu_x_sm{k+1}, u_meas_upd, t_upd, KF.param);
            h_u = h_u_(mu_x_sm{k+1}, u_meas_upd, t_upd, KF.param);
            h = h_(mu_x_sm{k+1}, u_meas_upd, t_upd, KF.param);
            
            R = Sigma2_z + h_u * Sigma2_u * h_u' + Sigma2_v;
            inv_R = inv(R);
            
            I_x_obs = h_x' * inv_R * h_x;
            i_x_obs = h_x' * inv_R * (z_meas_upd - (h - h_x * mu_x_sm{k+1}));
            
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
                u_meas_prev = u_meas{k-1};
                t_prev = t{k-1};
            else
                mu_x_nominal_prev = [];
                u_meas_prev = [];
                t_prev = [];
            end
            
            [mu_x_sm{k}, Sigma2_x_sm{k}, I_x_pred_bck, i_x_pred_bck] = IEKS_step(...
                k, mu_x_nominal_k, mu_x_nominal_prev, mu_x_fwd{k}, Sigma2_x_fwd{k}, ...
                u_meas{k}, z_meas{k}, t{k}, u_meas_prev, t_prev, I_x_pred_bck, i_x_pred_bck, ...
                KF.Delta_t, KF.Sigma_w, KF.Sigma_v, KF.Sigma_u, KF.Sigma_z, KF.param);
        end
    end
    
    % Store final iterated smoothed results in FilterResults data
    FilterResults.data.mu_x_sm = mu_x_sm;
    FilterResults.data.Sigma2_x_sm = Sigma2_x_sm;
    
    % Append names for unpack_simulation to dynamically pick up
    if ~any(strcmp(FilterResults.names, 'mu_x_sm'))
        FilterResults.names{end+1} = 'mu_x_sm';
    end
    if ~any(strcmp(FilterResults.names, 'Sigma2_x_sm'))
        FilterResults.names{end+1} = 'Sigma2_x_sm';
    end
end
