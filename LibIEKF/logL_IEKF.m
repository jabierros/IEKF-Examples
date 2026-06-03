function logL = logL_IEKF(S)
% LOGL_IEKF Computes the negative log-likelihood or joint energy function.
%
% logL = logL_IEKF(S)

    if ~isfield(S, 'smoother_in_likelihood')
        S.smoother_in_likelihood = 'none';
    end

    Sigma2_w_x = get_cov(S.sigma_w_x);
    Sigma2_u = get_cov(S.sigma_u);
    Sigma2_z = get_cov(S.sigma_z);
    Sigma2_v_x = get_cov(S.sigma_v_x);

    if strcmp(S.smoother_in_likelihood, 'none')
        % Initialize simulation state for standard forward run
        S.t = S.t_0;
        S.t_prev = S.t_0;
        S.x_true = S.x_true_0;
        S.x_true_prev = S.x_true_0;
        S.mu_x = S.mu_x_0;
        S.Sigma2_x = diag((S.sigma_x_0).^2);

        % Get initial measurements
        [u_meas, S] = get_u(S);
        S.u_meas = u_meas;
        [z_meas, S] = get_z(S);
        S.z_meas = z_meas;

        logL = 0;
        N_steps = S.t_end/S.Delta_t;
        
        for k = 1:N_steps
            S = IEKF_step(S);
            
            % Log-Likelihood computation based on forward prediction error
            h_x_ = h_x(S.mu_x_pred, S.u_meas, S.t, S.param);
            Sigma2_Innovation = h_x_ * S.Sigma2_x_pred * h_x_' + S.R;
            Innovation = S.z_meas - h(S.mu_x_pred, S.u_meas, S.t, S.param);
            
            logL = logL + 0.5 * (log(det(Sigma2_Innovation)) + Innovation' * inv(Sigma2_Innovation) * Innovation);
        end
    else
        % Run full forward pass with minimum logging required
        datalogging_string = {'t'; 'mu_x'; 'Sigma2_x'; 'u_meas'; 'z_meas'};
        Simulation = IEKF(S, datalogging_string);
        
        % Run chosen smoother
        if strcmp(S.smoother_in_likelihood, 'IEKS')
            Simulation = IEKS(Simulation, S);
        elseif strcmp(S.smoother_in_likelihood, 'IIEKS')
            Simulation = IIEKS(Simulation, S, 4);
        else
            error('Invalid value for S.smoother_in_likelihood. Must be ''none'', ''IEKS'', or ''IIEKS''.');
        end
        
        N = Simulation.k;
        t = Simulation.data.t;
        u_meas = Simulation.data.u_meas;
        z_meas = Simulation.data.z_meas;
        mu_x_sm = Simulation.data.mu_x_sm;
        
        logL = 0;
        
        for k = 1:N
            u_meas_curr = u_meas{k};
            z_meas_curr = z_meas{k};
            t_curr = t{k};
            
            % Measurement residual energy
            h_u_val = h_u(mu_x_sm{k}, u_meas_curr, t_curr, S.param);
            R = Sigma2_z + h_u_val * Sigma2_u * h_u_val' + Sigma2_v_x;
            Innovation = z_meas_curr - h(mu_x_sm{k}, u_meas_curr, t_curr, S.param);
            logL = logL + 0.5 * (Innovation' * inv(R) * Innovation + log(det(R)));
            
            % Transition energy from k-1 to k
            if k > 1
                u_meas_prev = u_meas{k};
                t_prev = t{k-1};
                f_u_val = f_u(mu_x_sm{k-1}, u_meas_prev, t_prev, S.param);
                Q = f_u_val * Sigma2_u * f_u_val' + Sigma2_w_x;
                Transition_Error = mu_x_sm{k} - f(mu_x_sm{k-1}, u_meas_prev, t_prev, S.param);
                logL = logL + 0.5 * (Transition_Error' * inv(Q) * Transition_Error + log(det(Q)));
            end
        end
    end
end

function cov_matrix = get_cov(sigma)
    if size(sigma, 2) == 1
        cov_matrix = diag(sigma.^2);
    else
        cov_matrix = sigma * sigma';
    end
end