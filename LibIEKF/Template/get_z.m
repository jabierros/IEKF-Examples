function [z_meas, S] = get_z(S)
% GET_Z Simulates sensor measurement with noise.
%
% [z_meas, S] = get_z(S)

    [x_true, S] = get_x_true(S);
    u_true = S.u_true_func(S.t);
    z_true = h_true(x_true, u_true, S.t, S.param_true);
    z_meas = z_true + S.sigma_z_true .* randn(size(z_true, 1), 1);
    
    S.x_true = x_true;
    S.u_true = u_true;
    S.z_true = z_true;
end