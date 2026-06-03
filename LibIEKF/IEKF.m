function Simulation = IEKF(S, datalogging_string)
% IEKF Performs the full forward pass of the Iterated Extended Kalman Filter.
%
% Simulation = IEKF(S, datalogging_string)

    N_steps = S.t_end/S.Delta_t;
    Simulation = init_simulation(N_steps + 1, datalogging_string, S);
    
    % Initialize simulation state
    S.t = S.t_0;
    S.t_prev = S.t_0;
    S.x_true = S.x_true_0;
    S.x_true_prev = S.x_true_0;
    S.mu_x = S.mu_x_0;
    S.Sigma2_x = S.Sigma2_x_0;
    
    [u_meas, S] = get_u(S); 
    S.u_meas = u_meas;
    [z_meas, S] = get_z(S);
    S.z_meas = z_meas;
    
    Simulation = log_simulation(Simulation, S);
    
    for k = 1:N_steps
        S = IEKF_step(S);
        Simulation = log_simulation(Simulation, S);
    end
end