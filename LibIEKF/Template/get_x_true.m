function [x_true, TrueSystem, SimOpts] = get_x_true(TrueSystem, SimOpts)
% GET_X_TRUE Integrates system equations using ode45.
%
% [x_true, TrueSystem, SimOpts] = get_x_true(TrueSystem, SimOpts)

    if SimOpts.t == SimOpts.t_0
        x_true = TrueSystem.x_true_0;
        TrueSystem.x_true_prev = x_true;
        SimOpts.t_prev = SimOpts.t_0;
    elseif SimOpts.t == SimOpts.t_prev
        x_true = TrueSystem.x_true_prev;
    elseif SimOpts.t > SimOpts.t_prev
        ode45_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
        [time_steps, x_steps] = ode45(@odefun, [SimOpts.t_prev, SimOpts.t], TrueSystem.x_true_prev, ode45_options);
        x_true = x_steps(end, :)';
        TrueSystem.x_true_prev = x_true;
        SimOpts.t_prev = SimOpts.t;
    else
        error('t must be such t==t_0, t=t_prev or t=t_prev+Delta_t');
    end

    function dxdt = odefun(t, x_)
        dxdt = dstate_true_(x_, u_true_func_(t), t, TrueSystem.param_true);
    end
end
