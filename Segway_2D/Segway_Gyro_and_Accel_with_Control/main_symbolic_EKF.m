clear all
close all

syms t x dx ddx theta dtheta ddtheta R delta1 delta3 m_Wh m_Ch I_Wh_G_Wh I_Ch_G_Ch c g Delta_t theta0 k_c c_c w_x_1 v_x_1 real

q=[x,theta]'
dq=[dx,dtheta]'
ddq=[ddx,ddtheta]'

OM_Wh=[0,dx/R,0]' %123 or xyz
AL_Wh=[0,ddx/R,0]' %123 or xyz
OM_Ch=[0,dtheta,0]' %123 or xyz
AL_Ch=[0,ddtheta,0]' %123 or xyz
Vel_C=[dx,0,0]' %xyz
Acc_C=[ddx,0,0]' %xyz
Vel_G_Wh=Vel_C
Acc_G_Wh=Acc_C
Vel_G_Ch=[dx*cos(theta),0,dx*sin(theta)]'+cross([0,dtheta,0]',[delta1,0,delta3]') %123
Acc_G_Ch=[ddx*cos(theta),0,+ddx*sin(theta)]'+cross(OM_Ch,cross(OM_Ch,[delta1,0,delta3]'))+ cross(AL_Ch,[delta1,0,delta3]') %123

iF_Wh=-m_Wh*Acc_C %xyz
iM_Wh_G_Wh=-I_Wh_G_Wh*AL_Wh %xyz
gF_Wh=[0,m_Wh*g*sin(theta),0] %xyz

iF_Ch=-m_Ch*Acc_G_Ch %123
iM_Ch_G_Ch=-I_Ch_G_Ch*AL_Ch %123
gF_Ch=[m_Ch*g*sin(theta),0,-m_Ch*g*cos(theta)] %123

M_m=k_c*(theta-theta0) +c_c*(dtheta);
M_M_P=[0,M_m,0]';%123 or xyz

M_vis=-c*OM_Wh

eq=(iF_Wh'*jacobian(Vel_G_Wh,dq))'+(iM_Wh_G_Wh'*jacobian(OM_Wh',dq))'+...
(gF_Wh*jacobian(Vel_G_Wh,dq))'+...
(iF_Ch'*jacobian(Vel_G_Ch,dq))'+(iM_Ch_G_Ch'*jacobian(OM_Ch,dq))'+...
(gF_Ch*jacobian(Vel_G_Ch,dq))'+...
(M_vis'*jacobian(OM_Wh',dq))'+...
(M_M_P'*jacobian(OM_Wh',dq))'-...
(M_M_P'*jacobian(OM_Ch',dq))'; eq=simplify(eq)

M=jacobian(eq,ddq); M=simplify(M)
delta=-subs(eq,ddq,[0,0]'); delta=simplify(delta)

ddq_aux=inv(M)*delta; ddq_aux=simplify(ddq_aux)

% State as a 1rst order diff eq.
x_=[q;dq]

% Input
u=[theta0]

% In this example we are not using noise source w_x_1, but we make this a
% general template able to export of funtions allowing for w_x.
% matlabFunction does not allow a funtion variable to be empty.
% So when using the funtions we use w_x_1=0 and diag_Sigma_w_x=[0]
w_x=[w_x_1];
v_x=[v_x_1];
dstate=[dq;ddq_aux]; dstate=simplify(dstate)
dstate_noisy=dstate; dstate_noisy=simplify(dstate_noisy)
dstate=subs(dstate,w_x,0); dstate=simplify(dstate)
dstate_x=jacobian(dstate);

% Param
param=[R delta1 delta3 m_Wh m_Ch I_Wh_G_Wh I_Ch_G_Ch c k_c c_c g Delta_t ]'

%% Euler discretization
f_noisy=x_+dstate_noisy*Delta_t; %x_k+1=f(x_k)
f_noisy=simplify(f_noisy);

%% Exponential discretization
% See eq.(43) in Exponential Integration Schemes in Multibody Dynamics Javier Ros, Xabier Iriarte, Aitor Plaza, Jorge Angeles (Multibody Dynamics 2011)
% f_noisy=x_ + Delta_t * inv(dstate_x*Delta_t)*(expm(dstate_x*Delta_t) - diag(ones(size(x_,1))))*dstate_noisy;

%%
f_w_x=jacobian(f_noisy,w_x); f_w_x=simplify(f_w_x)

f=subs(f_noisy,w_x,[0]); f=simplify(f)

f_x=jacobian(f,x_); f_x=simplify(f_x)

f_u=jacobian(f,u); f_u=simplify(f_u)

%Accelerometer at C in the Platform 
a=[ddx*cos(theta)-(+g*sin(theta))
   ddx*sin(theta)-(-g*cos(theta))] %123

a_ddq=jacobian(a,ddq)

a_ddq_0=subs(a,ddq,zeros(size(ddq)))

a=a_ddq*ddq_aux+a_ddq_0

h_noisy=[dtheta
   a ]; h_noisy=simplify(h_noisy)

h_v_x=jacobian(h_noisy,v_x); h_w_x=simplify(h_v_x)

h=subs(h_noisy,v_x,0); h=simplify(h)

h_x=jacobian(h,x_)

h_u=jacobian(h,u);h_u=simplify(h_u)

%%Export

matlabFunction(f,'file','f','vars',{x_ u t param});

matlabFunction(f_noisy,'file','f_noisy','vars',{x_ u t param w_x});

matlabFunction(f_x,'file','f_x','vars',{x_ u t param});

matlabFunction(f_u,'file','f_u','vars',{x_ u t param});

matlabFunction(f_w_x,'file','f_w_x','vars',{x_ u t param});

matlabFunction(h,'file','h','vars',{x_ u t param});

matlabFunction(h_noisy,'file','h_noisy','vars',{x_ u t param v_x});

matlabFunction(h_x,'file','h_x','vars',{x_ u t param});

matlabFunction(h_u,'file','h_u','vars',{x_ u t param});

matlabFunction(h_v_x,'file','h_v_x','vars',{x_ u t param});

matlabFunction(dstate,'file','dstate','vars',{x_ u t param});

matlabFunction(dstate_noisy,'file','dstate_noisy','vars',{x_ u t param w_x});

matlabFunction(dstate_x,'file','dstate_x','vars',{x_ u t param});

%% Disambiguation
% As the continuous equation is integrated in time to obtain the discrete 
% process eq. the continuous Wiener process w_x_1 is also integrated.
% Leaving appart the preceeding continuous equations
% in this script the meaning of the w_x_1 variable is a discrete
% Wiener proccess. This will be the meaning from now on, here and in the 
% numeric (discrete) counterpart of this listing. The relationship
% between the variance of the continuous Wiener procces and that of the
% discrete one is
% sqrt(Var(discrete_w_x_1))*Delta_t=sqrt(1/3*Delta_t^3*Var(continuous_w_x_1))=
% Delta_t*sqrt(1/3*Delta_t*Var((continuous_w_x_1)).
% Thus
%   sqrt(Var(discrete_w_x_1))=sqrt(1/3*Delta_t*Var((continuous_w_x_1)),
% and therefore
%   Var(discrete_w_x_1)=1/3*Delta_t*Var(continuous_w_x_1)
% But we don't need to be aware of this. As we only need to define
% properties for the discrete w_x_1 proccess.
% When integrating the continuous equation to generate the state series
% that will be used to generate the observation or measurement series
% Noise will be asumed to have a zero order hold. That is it will be
% constant during the integration step.