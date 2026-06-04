function FilterResults = IEKS(FilterResults, KF, SimOpts)
% IEKS Performs backward pass and smoothing fusion of the Information Extended Kalman Smoother.
%
% FilterResults = IEKS(FilterResults, KF, SimOpts)

    N = FilterResults.k;
    if N <= 0
        return;
    end
    
    n_x = size(KF.mu_x_0, 1);
    
    % Extract logged data
    t = FilterResults.data.t;
    u_meas = FilterResults.data.u_meas;
    z_meas = FilterResults.data.z_meas;
    mu_x_fwd = FilterResults.data.mu_x;
    Sigma2_x_fwd = FilterResults.data.Sigma2_x;
    
    % Preallocate cells for smoothed estimates
    mu_x_sm = cell(N, 1);
    Sigma2_x_sm = cell(N, 1);
    
    % Initialize backward pass variables
    I_x_pred_bck = zeros(n_x, n_x);
    i_x_pred_bck = zeros(n_x, 1);
    
    % Backward recursion loop
    for k = N:-1:1
        mu_x_nominal_k = mu_x_fwd{k};
        if k > 1
            mu_x_nominal_prev = mu_x_fwd{k-1};
            u_meas_prev = u_meas{k-1}; % input at k-1 used for k-1 -> k
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
    
    % Store smoothed results in FilterResults data
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
