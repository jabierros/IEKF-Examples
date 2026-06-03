function Simulation = log_simulation(Simulation, S)
% LOG_SIMULATION Logs state variables at the current time step.
%
% Simulation = log_simulation(Simulation, S)

    Simulation.k = Simulation.k + 1;
    k = Simulation.k;
    
    % Ensure sigma_x is computed for logging
    S.sigma_x = diag(S.Sigma2_x).^0.5;
    
    for i = 1:length(Simulation.names)
        name = Simulation.names{i};
        Simulation.data.(name){k} = S.(name);
    end
end
