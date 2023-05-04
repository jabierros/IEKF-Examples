function u_=get_u()
global t
global u_actual
global sigma_u_actual u_actual_func 
set_x_actual();
u_actual=u_actual_func(t);
% If there is an input it is supossed to be measured, measurement is
% generated from true value of u u_actual (obtained from the user provided
% function u_actual_func(t), to which sensor noise is added
u_=u_actual+sigma_u_actual.*randn(size(u_actual,1),1);
end