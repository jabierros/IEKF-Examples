function z_meas=get_z()
global t
global x_actual z_actual u_actual u_actual_func 
global Sigma_z_actual
global param
set_x_actual()
u_actual=u_actual_func(t);
%Sensor data is generated from true system estate and input (x,u), by using
%h, to which sensor noise is added
z_actual=h(x_actual,u_actual,t,param);
n_z=size(z_actual,1);
z_meas=z_actual+Sigma_z_actual*randn(n_z,1);
end