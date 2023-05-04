function z_=get_z()
global t
global x_actual z_actual u_actual u_actual_func 
global sigma_z_actual
global param
set_x_actual();
u_actual=u_actual_func(t);
% Sensor data is generated from true system state and input (x_actual,u_actual), by using
% h, to which sensor noise is added
z_actual=h(x_actual,u_actual,t,param);
z_=z_actual+sigma_z_actual.*randn(size(z_actual,1),1);
end