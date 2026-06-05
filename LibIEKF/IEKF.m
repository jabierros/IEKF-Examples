function FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string)
% IEKF Performs the full forward pass of the Iterated Extended Kalman Filter.
%
% FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string)

    N_steps = SimOpts.t_end / KF.Delta_t;
    FilterResults = init_simulation(N_steps + 1, datalogging_string, KF);
    
    % Initialize simulation state
    SimOpts.t = SimOpts.t_0;
    SimOpts.t_prev = SimOpts.t_0;
    TrueSystem.x_true = TrueSystem.x_true_0;
    TrueSystem.x_true_prev = TrueSystem.x_true_0;
    
    % Initialize filter state
    mu_x = KF.mu_x_0;
    if isvector(KF.sigma_x_0)
        Sigma2_x = diag(KF.sigma_x_0.^2);
    else
        Sigma2_x = KF.sigma_x_0 * KF.sigma_x_0';
    end
    
    % Initial epoch (t=0)
    [u_meas, TrueSystem, SimOpts] = meas_u(TrueSystem, SimOpts); 
    [z_meas, TrueSystem, SimOpts] = meas_z(TrueSystem, SimOpts);
    
    % Log initial epoch
    S_temp = struct();
    S_temp.t = SimOpts.t;
    S_temp.x_true = TrueSystem.x_true;
    S_temp.mu_x = mu_x;
    S_temp.Sigma2_x = Sigma2_x;
    if isfield(TrueSystem, 'u_true'), S_temp.u_true = TrueSystem.u_true; end
    S_temp.u_meas = u_meas;
    if isfield(TrueSystem, 'z_true'), S_temp.z_true = TrueSystem.z_true; end
    S_temp.z_meas = z_meas;
    S_temp.sigma_x = sqrt(diag(Sigma2_x));
    FilterResults = log_simulation(FilterResults, S_temp);
    
    for k = 1:N_steps
        % Keep previous input & time
        t_prev = SimOpts.t;
        u_meas_prev = u_meas;
        
        % Advance time
        SimOpts.t = t_prev + KF.Delta_t;
        SimOpts.t_prev = t_prev;
        
        % Read input and measurement at t (epoch k+1)
        [u_meas, TrueSystem, SimOpts] = meas_u(TrueSystem, SimOpts);
        [z_meas, TrueSystem, SimOpts] = meas_z(TrueSystem, SimOpts);
        
        % Run flat step
        [mu_x, Sigma2_x] = IEKF_step(...
            mu_x, Sigma2_x, u_meas_prev, u_meas, z_meas, t_prev, SimOpts.t, KF.Delta_t, ...
            KF.Sigma_w, KF.Sigma_v, KF.Sigma_u, KF.Sigma_z, KF.param);
            
        % Log current epoch
        S_temp = struct();
        S_temp.t = SimOpts.t;
        S_temp.x_true = TrueSystem.x_true;
        S_temp.mu_x = mu_x;
        S_temp.Sigma2_x = Sigma2_x;
        if isfield(TrueSystem, 'u_true'), S_temp.u_true = TrueSystem.u_true; end
        S_temp.u_meas = u_meas;
        if isfield(TrueSystem, 'z_true'), S_temp.z_true = TrueSystem.z_true; end
        S_temp.z_meas = z_meas;
        S_temp.sigma_x = sqrt(diag(Sigma2_x));
        FilterResults = log_simulation(FilterResults, S_temp);
    end
end