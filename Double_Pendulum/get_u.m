function u_meas=get_u()
global t
global u_actual n_u
global Sigma_u_actual u_actual_func 
set_x_actual();
u_actual=u_actual_func(t);
u_meas=u_actual+Sigma_u_actual*randn(n_u,1);
end