function [x_true, S] = get_x_true(S)
% GET_X_TRUE Integrates system equations using ode45.
%
% [x_true, S] = get_x_true(S)

    if S.t == S.t_0
        x_true = S.x_true_0;
        S.x_true_prev = x_true;
        S.t_prev = S.t_0;
    elseif S.t == S.t_prev
        x_true = S.x_true_prev;
    elseif S.t > S.t_prev
        ode45_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
        [time_steps, x_steps] = ode45(@odefun, [S.t_prev, S.t], S.x_true_prev, ode45_options);
        x_true = x_steps(end, :)';
        S.x_true_prev = x_true;
        S.t_prev = S.t;
    else
        error('t must be such t==t_0, t=t_prev or t=t_prev+Delta_t');
    end

    function dxdt = odefun(t, x_)
        dxdt = dstate_true(x_, S.u_true_func(t), t, S.param_true);
    end
end
