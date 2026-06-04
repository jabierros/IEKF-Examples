function [z_meas, TrueSystem, SimOpts] = get_z(TrueSystem, SimOpts)
% GET_Z Simulates sensor measurement with noise.
%
% [z_meas, TrueSystem, SimOpts] = get_z(TrueSystem, SimOpts)

    [x_true, TrueSystem, SimOpts] = get_x_true(TrueSystem, SimOpts);
    u_true = u_true_func_(SimOpts.t);
    z_true = h_true_(x_true, u_true, SimOpts.t, TrueSystem.param_true);
    z_meas = z_true + TrueSystem.sigma_z_true .* randn(size(z_true, 1), 1);
    
    TrueSystem.x_true = x_true;
    TrueSystem.u_true = u_true;
    TrueSystem.z_true = z_true;
end