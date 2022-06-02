function u_meas=get_u()
global t
global u_actual
global sigma_u_actual u_actual_func 
set_x_actual();
u_actual=u_actual_func(t);
u_meas=u_actual+sigma_u_actual.*randn(size(u_actual,1),1);
end