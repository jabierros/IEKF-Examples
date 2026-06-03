function [u_meas, S] = get_u(S)
% GET_U Simulates the measured input with noise.
%
% [u_meas, S] = get_u(S)

    u_true = S.u_true_func(S.t);
    u_meas = u_true + S.sigma_u_true .* randn(size(u_true, 1), 1);
    S.u_true = u_true;
end