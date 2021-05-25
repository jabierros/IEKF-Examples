clear all
close all

syms t x dx ddx m k x0 c Delta_t f_ext w_x_1 v_x_1 real

q=[x]'
dq=[dx]'
ddq=[ddx]'

param=[m k x0 c Delta_t]'

iF=-m*[ddx,0,0]' %xyz

sdF=[-k*(x-x0),0,0]'+[-c*dx,0,0]' %xyz

extF=[f_ext,0,0]' %xyz

V_G=[dx,0,0]' %xyz

V_A=V_G %xyz

V_B=V_G %xyz

eq=-((iF'*jacobian(V_G,dq))'+...
   (sdF'*jacobian(V_A,dq))'+...
   (extF'*jacobian(V_B,dq))')

eq=simplify(eq)

M=jacobian(eq,ddq)
M=simplify(M)
delta=-subs(eq,ddq,[0]')

delta=simplify(delta)

ddq_aux=inv(M)*delta
ddq_aux=simplify(ddq_aux)

x_=[q;dq]
u=sym([f_ext])

% In this example we are not using noise source w_x_1, but we need it to
% make a general export of funtions that allows for w_x, matlabFunction
% does not allows funtion variable to be empty.
w_x=[w_x_1];
v_x=[v_x_1];

dstate=[dq;ddq_aux]; dstate=simplify(dstate)

dstate_noisy=dstate; dstate_noisy=simplify(dstate_noisy)
dstate=subs(dstate,w_x,[0]); dstate=simplify(dstate)
dstate_x=jacobian(dstate);

f=x_+dstate_noisy*Delta_t; f=simplify(f) %x_k+1=f(x_k)

f_noisy=f; % no special noise in this example, but using a programing template that assumes that it can exist

f_w_x=jacobian(f_noisy,w_x); f_w_x=simplify(f_w_x)

f=subs(f,w_x,[0]); f=simplify(f)

f_x=jacobian(f,x_); f_x=simplify(f_x)

f_u=jacobian(f,u); f_u=simplify(f_u)

%Position sensors in x and y axes
h=[ddq_aux]; h=simplify(h)

h_noisy=h; % no special noise in this example, but using a programing template that assumes that it can exist

h_v_x=jacobian(h_noisy,v_x); h_v_x=simplify(h_v_x)

h=subs(h,v_x,[0]); h=simplify(h)

h_x=jacobian(h,x_); h_x=simplify(h_x)

h_u=jacobian(h,u); h_u=simplify(h_u)

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
