function logL = logL_IEKF(KF, TrueSystem, SimOpts)
% LOGL_IEKF Computes the negative log-likelihood or joint energy function.
%
% logL = logL_IEKF(KF, TrueSystem, SimOpts)

    if ~isfield(KF, 'smoother_in_likelihood')
        KF.smoother_in_likelihood = 'none';
    end

    Sigma2_w = KF.Sigma_w * KF.Sigma_w';
    Sigma2_u = KF.Sigma_u * KF.Sigma_u';
    Sigma2_z = KF.Sigma_z * KF.Sigma_z';
    Sigma2_v = KF.Sigma_v * KF.Sigma_v';

    if strcmp(KF.smoother_in_likelihood, 'none')
        % Initialize simulation state for standard forward run
        SimOpts.t = SimOpts.t_0;
        SimOpts.t_prev = SimOpts.t_0;
        TrueSystem.x_true = TrueSystem.x_true_0;
        TrueSystem.x_true_prev = TrueSystem.x_true_0;
        mu_x = KF.mu_x_0;
        if isvector(KF.sigma_x_0)
            Sigma2_x = diag(KF.sigma_x_0.^2);
        else
            Sigma2_x = KF.sigma_x_0 * KF.sigma_x_0';
        end

        % Get initial measurements
        [u_meas, TrueSystem, SimOpts] = get_u(TrueSystem, SimOpts);
        [z_meas, TrueSystem, SimOpts] = get_z(TrueSystem, SimOpts);

        logL = 0;
        N_steps = SimOpts.t_end/KF.Delta_t;
        
        for k = 1:N_steps
            t_prev = SimOpts.t;
            u_meas_prev = u_meas;
            
            SimOpts.t = t_prev + KF.Delta_t;
            SimOpts.t_prev = t_prev;

            [u_meas, TrueSystem, SimOpts] = get_u(TrueSystem, SimOpts);
            [z_meas, TrueSystem, SimOpts] = get_z(TrueSystem, SimOpts);

            % Call step filter to propagate and update
            [mu_x_next, Sigma2_x_next, mu_x_pred, Sigma2_x_pred] = IEKF_step(...
                mu_x, Sigma2_x, u_meas_prev, u_meas, z_meas, t_prev, SimOpts.t, KF.Delta_t, ...
                KF.Sigma_w, KF.Sigma_v, KF.Sigma_u, KF.Sigma_z, KF.param);
            
            % Log-Likelihood computation based on forward prediction error
            h_x = h_x_(mu_x_pred, u_meas, SimOpts.t, KF.param);
            h_u = h_u_(mu_x_pred, u_meas, SimOpts.t, KF.param);
            R = Sigma2_z + h_u * Sigma2_u * h_u' + Sigma2_v;
            
            Sigma2_Innovation = h_x * Sigma2_x_pred * h_x' + R;
            Innovation = z_meas - h_(mu_x_pred, u_meas, SimOpts.t, KF.param);
            
            logL = logL + 0.5 * (log(det(Sigma2_Innovation)) + Innovation' * inv(Sigma2_Innovation) * Innovation);
            
            mu_x = mu_x_next;
            Sigma2_x = Sigma2_x_next;
        end
    else
        % Run full forward pass with minimum logging required
        datalogging_string = {'t'; 'mu_x'; 'Sigma2_x'; 'u_meas'; 'z_meas'};
        FilterResults = IEKF(KF, TrueSystem, SimOpts, datalogging_string);
        
        % Run chosen smoother
        if strcmp(KF.smoother_in_likelihood, 'IEKS')
            FilterResults = IEKS(FilterResults, KF, SimOpts);
        elseif strcmp(KF.smoother_in_likelihood, 'IIEKS')
            FilterResults = IIEKS(FilterResults, KF, SimOpts, 4);
        else
            error('Invalid value for KF.smoother_in_likelihood. Must be ''none'', ''IEKS'', or ''IIEKS''.');
        end
        
        N = FilterResults.k;
        t = FilterResults.data.t;
        u_meas = FilterResults.data.u_meas;
        z_meas = FilterResults.data.z_meas;
        mu_x_sm = FilterResults.data.mu_x_sm;
        
        logL = 0;
        
        for k = 1:N
            u_meas_curr = u_meas{k};
            z_meas_curr = z_meas{k};
            t_curr = t{k};
            
            % Measurement residual energy
            h_u = h_u_(mu_x_sm{k}, u_meas_curr, t_curr, KF.param);
            R = Sigma2_z + h_u * Sigma2_u * h_u' + Sigma2_v;
            Innovation = z_meas_curr - h_(mu_x_sm{k}, u_meas_curr, t_curr, KF.param);
            logL = logL + 0.5 * (Innovation' * inv(R) * Innovation + log(det(R)));
            
            % Transition energy from k-1 to k
            if k > 1
                u_meas_prev = u_meas{k-1};
                t_prev = t{k-1};
                f_u = f_u_(mu_x_sm{k-1}, u_meas_prev, t_prev, KF.param, KF.Delta_t);
                Q = f_u * Sigma2_u * f_u' + Sigma2_w;
                Transition_Error = mu_x_sm{k} - f_(mu_x_sm{k-1}, u_meas_prev, t_prev, KF.param, KF.Delta_t);
                logL = logL + 0.5 * (Transition_Error' * inv(Q) * Transition_Error + log(det(Q)));
            end
        end
    end
end