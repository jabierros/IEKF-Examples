function unpack_simulation(Simulation)
% UNPACK_SIMULATION Unpacks logged variables into caller's workspace as series.
%
% unpack_simulation(Simulation)

    N_logged = Simulation.k;
    
    for i = 1:length(Simulation.names)
        name = Simulation.names{i};
        cell_data = Simulation.data.(name)(1:N_logged);
        if N_logged > 0
            sz = size(Simulation.data.(name){1});
        else
            sz = [0, 0];
        end
        
        if prod(sz) == 0
            % Handle empty matrices (e.g., when n_u = 0)
            series_data = zeros(N_logged, 0);
        elseif sz(1) == 1 || sz(2) == 1
            % Scalar or Vector: stack into a 2D matrix of size N_logged x max(sz)
            % to be backward-compatible with 2D series (time along rows)
            series_data = zeros(N_logged, max(sz));
            for k = 1:N_logged
                series_data(k, :) = cell_data{k}(:).';
            end
        else
            % Matrix (n x m): stack into a 3D matrix of size n x m x N_logged
            series_data = zeros(sz(1), sz(2), N_logged);
            for k = 1:N_logged
                series_data(:, :, k) = cell_data{k};
            end
        end
        
        assignin('caller', [name, '_series'], series_data);
    end
end
