function [u_meas, TrueSystem, SimOpts] = get_u(TrueSystem, SimOpts)
% GET_U Simulates the measured input with noise.
%
% [u_meas, TrueSystem, SimOpts] = get_u(TrueSystem, SimOpts)

    u_true = u_true_func_(SimOpts.t);
    u_meas = u_true + TrueSystem.sigma_u_true .* randn(size(u_true, 1), 1);
    TrueSystem.u_true = u_true;
end