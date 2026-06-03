function Simulation = init_simulation(N, datalogging_string, S)
% INIT_SIMULATION Preallocates cell arrays for in-memory logging.
%
% Simulation = init_simulation(N, datalogging_string, S)

    Simulation = struct();
    Simulation.N = N;
    Simulation.k = 0;
    Simulation.names = datalogging_string;
    Simulation.data = struct();
    
    for i = 1:length(datalogging_string)
        name = datalogging_string{i};
        Simulation.data.(name) = cell(N, 1);
    end
end
