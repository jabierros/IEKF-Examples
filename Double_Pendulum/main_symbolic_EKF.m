clear all; clc;
% cd(fileparts(which(mfilename)));  % In a .m file uncomment this
cd(fileparts(matlab.desktop.editor.getActiveFilename)); % In a .mlx file uncomment this
%% Example: Double Pendulum

syms t...
     theta1 dtheta1 ddtheta1...
     theta2 dtheta2 ddtheta2...
     m1 m2 l1 l2 g Delta_t     real
param=[m1 m2 l1 l2 g]';
param_true=param;
%% 

q=[theta1, theta2]';
dq=[dtheta1, dtheta2]';
ddq=[ddtheta1, ddtheta2]';
%% 


acc_P1=[cos(theta1), 0,-sin(theta1);
        0,           1, 0;
        sin(theta1), 0, cos(theta1)]*[ddtheta1*l1;
                                      0
                                      dtheta1^2*l1];%xyz
    
acc_P2=acc_P1+[cos(theta2), 0,-sin(theta2);
           0,           1, 0;
           sin(theta2), 0, cos(theta2)]*[ddtheta2*l2;
                                         0
                                         dtheta2^2*l2];%xyz
                               
iF_P1=-m1*acc_P1; %xyz
iF_P2=-m2*acc_P2; %xyz
gF_P1=[0,0,-m1*g]'; %xyz
gF_P2=[0,0,-m2*g]'; %xyz

v_P1=[cos(theta1), 0,-sin(theta1);
      0,           1, 0;
      sin(theta1), 0, cos(theta1)]*[dtheta1*l1;
                                    0
                                    0]; %xyz
v_P2=v_P1+[cos(theta2), 0,-sin(theta2);
           0,           1, 0;
           sin(theta2), 0, cos(theta2)]*[dtheta2*l2;
                                         0
                                         0]; %xyz
                               
Dyn_eq=-((iF_P1'*jacobian(v_P1,dq))'+...
     (iF_P2'*jacobian(v_P2,dq))'+...
     (gF_P1'*jacobian(v_P1,dq))'+...
     (gF_P2'*jacobian(v_P2,dq))'); Dyn_eq=simplify(Dyn_eq)
%% 
% Virtual power dynamic equations structure

M_qq=jacobian(Dyn_eq,ddq); M_qq=simplify(M_qq)
delta_q=-subs(Dyn_eq,ddq,zeros(size(ddq))); delta_q=simplify(delta_q)
%% 

ddq_func=inv(M_qq)*delta_q; ddq_func=simplify(ddq_func)
x_=[q;dq]
%% 

u_=sym(zeros(0,1))
u_true=sym(zeros(0,1))
%% 

dstate=[dq;ddq_func]; dstate=simplify(dstate)
%% 

dstate_x=jacobian(dstate,x_)
f=x_+dstate*Delta_t; f=simplify(f) %x_k+1=f(x_k)
f_x=jacobian(f,x_); f_x=simplify(f_x)
f_u=jacobian(f,u_); f_u=simplify(f_u)
%% 

acc_P2_Btheta2=simplify(subs([cos(theta2), 0,-sin(theta2);
                              0,           1, 0;
                              sin(theta2), 0, cos(theta2)]'*(acc_P2-[0,0,-g]'),...
                             ddq,ddq_func));
h=[dtheta2;
   acc_P2_Btheta2'*[1,0,0]'
   acc_P2_Btheta2'*[0,0,1]']; h=simplify(h)
h_x=jacobian(h,x_); h_x=simplify(h_x)
h_u=jacobian(h,u_); h_u=simplify(h_u)
h_true=h;

x_true=x_;
u_true=u_;
dstate_true=dstate;

matlabFunction(f,'file','f_','vars',{x_, u_, t, param, Delta_t});
matlabFunction(f_x,'file','f_x_','vars',{x_, u_, t, param, Delta_t});
matlabFunction(f_u,'file','f_u_','vars',{x_, u_, t, param, Delta_t});
matlabFunction(h,'file','h_','vars',{x_, u_, t, param});
matlabFunction(h_x,'file','h_x_','vars',{x_, u_, t, param});
matlabFunction(h_u,'file','h_u_','vars',{x_, u_, t, param});
matlabFunction(h_true,'file','h_true_','vars',{x_true, u_true, t, param_true});
matlabFunction(dstate,'file','dstate_','vars',{x_, u_, t, param});
matlabFunction(dstate_x,'file','dstate_x_','vars',{x_, u_, t, param});
matlabFunction(dstate_true,'file','dstate_true_','vars',{x_true, u_true, t, param_true});
matlabFunction(u_true,'file','u_true_func_','vars',{t});