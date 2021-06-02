clear all
close all

syms t...
    theta1 dtheta1 ddtheta1...
    theta2 dtheta2 ddtheta2...
    m1 m2 l1 l2 g Delta_t...
    w_x_1 v_x_1 real

q=[theta1, theta2]'
dq=[dtheta1, dtheta2]'
ddq=[ddtheta1, ddtheta2]'

param=[m1 m2 l1 l2 g Delta_t]'

Acc_P1=[cos(theta1), 0,-sin(theta1);
        0,           1, 0;
        sin(theta1), 0, cos(theta1)]*[ddtheta1*l1;
                                      0
                                      dtheta1^2*l1];%xyz
    
Acc_P2=Acc_P1+[cos(theta2), 0,-sin(theta2);
           0,           1, 0;
           sin(theta2), 0, cos(theta2)]*[ddtheta2*l2;
                                         0
                                         dtheta2^2*l2];%xyz
                               
iF_P1=-m1*Acc_P1 %xyz
iF_P2=-m2*Acc_P2 %xyz

V_P1=[cos(theta1), 0,-sin(theta1);
      0,           1, 0;
      sin(theta1), 0, cos(theta1)]*[dtheta1*l1;
                                    0
                                    0]; %xyz
V_P2=V_P1+[cos(theta2), 0,-sin(theta2);
           0,           1, 0;
           sin(theta2), 0, cos(theta2)]*[dtheta2*l2;
                                         0
                                         0]; %xyz
                               
gF_P1=[0,0,-m1*g]';%xyz
gF_P2=[0,0,-m2*g]';%xyz
eq=-((iF_P1'*jacobian(V_P1,dq))'+...
     (iF_P2'*jacobian(V_P2,dq))'+...
     (gF_P1'*jacobian(V_P1,dq))'+...
     (gF_P2'*jacobian(V_P2,dq))')

eq=simplify(eq)

M=jacobian(eq,ddq)
M=simplify(M)
delta=-subs(eq,ddq,zeros(size(ddq)))

delta=simplify(delta)

ddq_aux=inv(M)*delta
ddq_aux=simplify(ddq_aux)

x_=[q;dq]
u=sym([])

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
Acc_P2_Btheta2=simplify(subs([cos(theta2), 0,-sin(theta2);
                              0,           1, 0;
                              sin(theta2), 0, cos(theta2)]'*(Acc_P2-[0,0,-g]'),...
                             ddq,ddq_aux));
    
h=[dtheta2;
   Acc_P2_Btheta2'*[1,0,0]'
   Acc_P2_Btheta2'*[0,0,1]']; h=simplify(h)

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
