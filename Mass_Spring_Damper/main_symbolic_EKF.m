clear all; clc;
try
    cd(fileparts(matlab.desktop.editor.getActiveFilename));
catch
    try
        cd(fileparts(which(mfilename)));
    catch
        % Keep current directory
    end
end
%% Example: Mass Spring Damper

%% Solution:
%% This is the first part of the solution: Determine model and sensor equations, and export the required functions for simulation and Kalman filtering.

%% Parameters

syms t x dx ddx m k rho0 c Delta_t f_ext real
param=[m k rho0 c]';
param_true=param;
%% 

q=[x]';
dq=[dx]';
ddq=[ddx]';
%% 

iF=-m*[ddx,0,0]' %xyz
sdF=[-k*(x-rho0),0,0]'+[-c*dx,0,0]' %xyz
extF=[f_ext,0,0]' %xyz

V_G=[dx,0,0]' %xyz
V_A=V_G %xyz
V_B=V_G %xyz

Dyn_eq=-((iF'*jacobian(V_G,dq))'+...
   (sdF'*jacobian(V_A,dq))'+...
   (extF'*jacobian(V_B,dq))'); Dyn_eq=simplify(Dyn_eq)
%% 

M_qq=jacobian(Dyn_eq,ddq); M_qq=simplify(M_qq)
delta_q=-subs(Dyn_eq,ddq,sym(0)*ddq); delta_q=simplify(delta_q)
%% 

ddq_func=inv(M_qq)*delta_q; ddq_func=simplify(ddq_func)
%% 

x_=[q;dq]
%% 

u=sym([f_ext])
%% 

dstate=[dq;ddq_func]; dstate=simplify(dstate)
%% 

dstate_x=jacobian(dstate,x_)
%% 

f=x_+dstate*Delta_t; f=simplify(f) %x_k+1=f(x_k)

f_x=jacobian(f,x_); f_x=simplify(f_x)

f_u=jacobian(f,u); f_u=simplify(f_u)

%% 

h=[ddq_func]; h=simplify(h)

h_x=jacobian(h,x_); h_x=simplify(h_x)

h_u=jacobian(h,u); h_u=simplify(h_u)

x_true=x_;
syms u_true real
dstate_true=subs(dstate,u,u_true);
h_true=subs(h,u,u_true);

u_true_func=sym(0);

matlabFunction(f,'file','f_','vars',{x_, u, t, param, Delta_t});
matlabFunction(f_x,'file','f_x_','vars',{x_, u, t, param, Delta_t});
matlabFunction(f_u,'file','f_u_','vars',{x_, u, t, param, Delta_t});
matlabFunction(h,'file','h_','vars',{x_, u, t, param});
matlabFunction(h_x,'file','h_x_','vars',{x_, u, t, param});
matlabFunction(h_u,'file','h_u_','vars',{x_, u, t, param});
matlabFunction(h_true,'file','h_true_','vars',{x_true, u_true, t, param_true});
matlabFunction(dstate,'file','dstate_','vars',{x_, u, t, param});
matlabFunction(dstate_x,'file','dstate_x_','vars',{x_, u, t, param});
matlabFunction(dstate_true,'file','dstate_true_','vars',{x_true, u_true, t, param_true});
matlabFunction(u_true_func,'file','u_true_func_','vars',{t});