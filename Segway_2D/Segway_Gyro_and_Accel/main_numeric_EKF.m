
clear all
close all
seed=1789
rng(seed) %every run random numbers are repeated
ode45_options=odeset('RelTol',1e-12); %,'AbsTol',1e-12)

% Segway Gyro and Accel

%% Parameter and initial state definition
R=0.25
delta1=0.0
delta3=1
m_R=15
m_P=80
I_R=1/2*(0.5*m_R)*R^2
I_P=m_P*(delta1^2+delta3^2)/(4)^2
c=20/10*0.25
g=9.8
Delta_t=0.001

fig_dir=['EKF_',num2str(Delta_t)];
mkdir(fig_dir);

param=[R delta1 delta3 m_R m_P I_R I_P c g Delta_t ]'

x=0
theta=0.01*pi/180
dx=0
dtheta=0

q=[x,theta]'
dq=[dx,dtheta]'
M_m=0

x_0_=[q;dq]
u_0_=[M_m]

x_=x_0_
u_=u_0_

T=30

%% Generating input and noise to run and validate the filter

t_series=[0:Delta_t:T]

u_series=25*(heaviside(t_series-5)-heaviside(t_series-15));

du_series=0*t_series

%Sensor white noise definition and generation
sigma_gyro=0.5; % rad/s
sigma_acc_x=0.1*10; % m/s^2
sigma_acc_z=0.1*10; % m/s^2

Sigma_z=diag([sigma_gyro,sigma_acc_x,sigma_acc_z])

n_z=size(Sigma_z,1)
z_series_noise=[];
for k=1:length(t_series)
    z_series_noise=[z_series_noise, Sigma_z*randn(n_z,1)];
end

%Input white noise definition and generation
sigma_u=1*0.1; % N m

n_u=size(u_,1)
u_series_noise=[];
for k=1:length(t_series)
    u_series_noise=[u_series_noise,  sigma_u*randn(n_u,1)];
end

%% Initial filter state

delta_x0=10 % m
delta_theta0=30*pi/180 % rad 
delta_dx0=1 %  m/s
delta_dtheta0=1/(delta1^2+delta3^2)^0.5 %angular velocity if at  1m/s the whhel stops abruptly

%initial state expected value
mu_x_0_=x_0_+[delta_x0;
       delta_theta0;
       delta_dx0;
       delta_dtheta0];

mu_x_=mu_x_0_;

%initial state covariance
Sigma2_x_0_ = diag([delta_x0^2,delta_theta0^2,delta_dx0^2,delta_dtheta0^2])
Sigma2_x_ = Sigma2_x_0_;



%% Read input (k=1)
t=t_series(1,1);
u_=u_series(:,1)+u_series_noise(:,1);
z_without_noise=h(x_,u_series(:,1),t,param);
z_=z_without_noise+z_series_noise(:,1);

% Effect of noise in u on model covariance (f_u_ is dependent on the state)
f_u_=f_u(mu_x_,u_,t,param);
Q_u= f_u_*sigma_u^2*f_u_'

% Discretization error effect on model covariance
Q_discr= (1/factorial(2)*diag([100*Delta_t^2,100*Delta_t^2,6000*Delta_t^2,6000*Delta_t^2]))^2;

% Discretization error and measured input error are assumed independent
Q=Q_discr+Q_u

h_u_=h_u(f(mu_x_,u_,t,param),u_series(1+1),t+Delta_t,param)
% Effect of noise in u on sensor equation covariance
R_u= h_u_*sigma_u^2*h_u_'

R_z=diag([sigma_gyro^2,sigma_acc_x^2,sigma_acc_z^2])
% Sensor error and measured input error are assumed independent
R=R_z+R_u

x_series=[x_];
mu_x_series=[mu_x_];
z_series=[h(x_,u_,t,param)];
z_series_noisy=[z_];
sigma_x_=diag(Sigma2_x_).^0.5;
sigma_x_series=[sigma_x_];
pred_error_series=[zeros(size(x_))];
K_series=[reshape(zeros(4,3),12,1)];

for k=1:T/Delta_t
    
%% Read input (k)
    t_prev=t;
    
%%---Begin--- One step of Kalman filter (k-th)

%% Prediction
 f_u_=f_u(mu_x_,u_,t,param);
% effect of noise in u on model covariance (f_u_ is dependent on the state)
 Q_u= f_u_*sigma_u^2*f_u_';

% Discretization error effect on model covariance
 Q_discr= (1/factorial(2)*diag([100*Delta_t^2,100*Delta_t^2,6000*Delta_t^2,6000*Delta_t^2]))^2;

% Discretization error and measured input error are assumed independent
 Q=Q_discr+Q_u;

F_=F(mu_x_,u_,t,param);
mu_x_pred_ = f(mu_x_,u_,t,param);
 Sigma2_x_pred = F_ * Sigma2_x_ * F_' + Q;

%% Read Input and Sensor (k+1)
 %Read input k+1
 u_=u_series(:,k+1)+u_series_noise(:,k+1);
 
%Read sensors k+1
    %z_=z_series(:,k+1)+z_series_noise(:,k+1);
    t_prev=t;

    %Sensor data is generated from true system estate (using h) and adding noise
    %True system state is generated by integration using MATLAB ode45 integrator 
t=t_series(1,k+1);
    [time_steps,x_steps] = ode45(@(t,x_) dstate(x_, u_series(:,k), t,param),[t_series(1,k),t_series(1,k+1)],x_,ode45_options);
    x_=x_steps(end,:)';
    z_without_noise=h(x_,u_series(:,k+1),t,param);
    z_=z_without_noise+z_series_noise(:,k+1);
    pred_error=f(x_series(:,k),u_series(:,k),t_series(1,k),param)-x_;
    %Exponential discretization A and B are continuous state Jacobian wrt 
    %x and u
    %x_=x_ + inv(A)*(expm(A*Delta_t) - diag(ones(size(x_,1))))*B u_;
 
%% Correction and Fussion
 h_u_= h_u(mu_x_pred_,u_,t,param) ;
% effect of noise in u on sensor equation covariance 
 R_u= h_u_*sigma_u^2*h_u_';
 
 R_z=diag([sigma_gyro^2,sigma_acc_x^2,sigma_acc_z^2]);
% Sensor error and measured input error are assumed independent
 R=R_z+R_u;

H_= H(mu_x_pred_,u_,t,param);
S = H_ * Sigma2_x_pred * H_' + R; % Compute Innovation Covariance
K_gain = Sigma2_x_pred * H_' * inv(S); % Compute  and Kalman Gain

mu_x_ = mu_x_pred_ + K_gain * ( z_ - h(mu_x_pred_,u_,t,param) );

%Sigma2_x_ = Sigma2_x_pred - K_gain * S * K_gain';
Sigma2_x_ = Sigma2_x_pred - K_gain * H_ * Sigma2_x_pred ; %A simplification of the previous from
%Sigma2_x_ = (diag(length(x_est_)) - K_gain * H_) * Sigma2_x_pred * (diag(length(x_est_)) - K_gain * H_)' + K_gain * R * K_gain'; %Joseph form

%%---End--- One step of Kalman filter (k-th)

    x_series=[x_series,x_];
mu_x_series=[mu_x_series, mu_x_];
    z_series=[z_series, z_without_noise];
    z_series_noisy=[z_series_noisy, z_];
sigma_x_series=[sigma_x_series, sqrt(diag(Sigma2_x_))];
pred_error_series=[pred_error_series, pred_error];
K_series=[K_series,reshape(K_gain,12,1)];

end

figure
plot(t_series,x_series)

I =legend('$x$','$\theta$','$\dot{x}$','$\dot{\theta}$')
set(I,'interpreter','latex');
I =title('Real state (not known in a real experiment)')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','trajectory'],'epsc')
saveas(gcf,[fig_dir,'/','trajectory'],'png')

figure
plot(t_series,u_series,'-',t_series,u_series+u_series_noise,'--');
I =legend('$\mathbf{u}=[M_m]$','${meas}(\mathbf{u})=[{meas}(M_m)]$')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','u'],'epsc')
saveas(gcf,[fig_dir,'/','u'],'png')

x_=x_0_
x_series_Euler=[]
for k=1:length(t_series)
    x_series_Euler=[x_series_Euler,x_];
    x_=f(x_, u_series(k), t_series(k),param);
end

figure
plot(t_series,x_series_Euler)
I =legend('$x$','$\theta$','$\dot{x}$','$\dot{\theta}$')
set(I,'interpreter','latex');
I =title('"Real state" as obtained by Euler discretization (not known in a real experiment)')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','trajectory_Euler'],'epsc')
saveas(gcf,[fig_dir,'/','trajectory_Euler'],'png')

figure
ddx_series=[];
for i=1:length(t_series)
    ddx_series=[ddx_series,ddstate(x_series(:,i),u_series(:,i),du_series(:,i),t_series(1,i),param)];
end

plot(t_series,ddx_series)
I =legend('$\ddot{x}$','$\ddot{\theta}$','$\dot{\ddot{x}}$','$\dot{\ddot{\theta}}$')
set(I,'interpreter','latex');
I =title('Approximation of $\ddot{\mathbf{x}}$ with $\dot{\mathbf{u}}=\mathbf{0}$ to estimate Euler $2^{nd}$ order discretization error')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','discrete_2nd_x_derivative'],'epsc')
saveas(gcf,[fig_dir,'/','discrete_2nd_x_derivative'],'png')

figure
plot(t_series,z_series,'-',t_series,z_series_noisy,'--')
I =legend('$\dot{\theta}$','$a_x$','$a_y$','${meas}(\dot{\theta})$','${meas}(a_x)$','${meas}(a_y)$')
set(I,'interpreter','latex');
I=title('$\mathbf{z}$ with and without noise')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','z_without_z_with_noise'],'epsc')
saveas(gcf,[fig_dir,'/','z_without_z_with_noise'],'png')

figure
plot(t_series,z_series_noise,'-')
I =legend('$\dot{\theta}-{meas}(\dot{\theta})$','$a_x-{meas}(a_x)$','$a_y-{meas}(a_y)$')
set(I,'interpreter','latex');
I=title('$\mathbf{z}-{meas}(\mathbf{z})$ with and without noise')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','z_with_noise_minus_z_without'],'epsc')
saveas(gcf,[fig_dir,'/','z_with_noise_minus_z_without'],'png')

figure
plot(t_series,x_series,'-',t_series,mu_x_series,'--')
xlim([0 T])
I = legend('$x$','$\theta$','$\dot{x}$','$\dot{\theta}$','$\hat{\mu}_{x}$','$\hat{\mu}_\theta$','$\hat{\mu}_{\dot{x}}$','$\hat{\mu}_{\dot{\theta}}$');
set(I,'interpreter','latex');
I=title('$\mathbf{x}$ vs $\hat{\mu}_{\mathbf{x}}$')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','mu_actual'],'epsc')
saveas(gcf,[fig_dir,'/','mu_actual'],'png')

figure
plot(t_series,x_series-mu_x_series,'-')
xlim([0 T])
I = legend('$x-\hat{\mu}_{x}$','$\theta-\hat{\mu}_\theta$','$\dot{x}-\hat{\mu}_{\dot{x}}$','$\dot{\theta}-\hat{\mu}_{\dot{\theta}}$');
set(I,'interpreter','latex');
I=title('$\hat{\mu}_{\mathbf{x}}-\mathbf{x}$')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','Error_mu_minus_actual'],'epsc')
saveas(gcf,[fig_dir,'/','Error_mu_minus_actual'],'png')

figure
semilogy(t_series,sigma_x_series,'-')
I =legend('$\sigma_{x}$','$\sigma_{\theta}$','$\sigma_{\dot{x}}$','$\sigma_{\dot{\theta}}$')
set(I,'interpreter','latex');
I=title('$diag({{\Sigma}}_{\mathbf{x}})$')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','sigma_x'],'epsc')
saveas(gcf,[fig_dir,'/','sigma_x'],'png')

figure
semilogy(t_series,abs(pred_error_series),'-')
I =legend('$\varepsilon_{x}$','$\varepsilon_{\theta}$','$\varepsilon_{\dot{x}}$','$\varepsilon_{\dot{\theta}}$')
set(I,'interpreter','latex');
I=title('Euler method prediction error')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','varepsilon_x'],'epsc')
saveas(gcf,[fig_dir,'/','varepsilon_x'],'png')

figure
plot(t_series,K_series,'-')
I =legend('$k_{11}$','$k_{12}$','$k_{13}$','$k_{21}$','$k_{22}$','$k_{23}$','$k_{31}$','$k_{32}$','$k_{33}$','$k_{41}$','$k_{42}$','$k_{43}$')
set(I,'interpreter','latex');
I=title('Kalman gain $\mathbf{K}$')
set(I,'interpreter','latex');
saveas(gcf,[fig_dir,'/','Kalman_gain'],'epsc')
saveas(gcf,[fig_dir,'/','Kalman_gain'],'png')

% Observability

%Linear observability near the end of the simulation
OB=obsv(F_,H_)

rank(OB)
size(OB) %Nonlinear observabilty matrix mus be bigger (more rows) than the linear one.

%Nonlinear observability near the end of the simulation
k=k-5

OB=[H(x_series(:,k),u_series(:,k),t_series(:,k),param),
    H(x_series(:,k+1),u_series(:,k+1),t_series(:,k+1),param)*F(x_series(:,k+1),u_series(:,k+1),t_series(:,k+1),param)
    H(x_series(:,k+2),u_series(:,k+2),t_series(:,k+2),param)*F(x_series(:,k+2),u_series(:,k+2),t_series(:,k+2),param)*F(x_series(:,k+1),u_series(:,k+1),t_series(:,k+1),param)
    H(x_series(:,k+3),u_series(:,k+3),t_series(:,k+3),param)*F(x_series(:,k+3),u_series(:,k+3),t_series(:,k+3),param)*F(x_series(:,k+2),u_series(:,k+2),t_series(:,k+2),param)*F(x_series(:,k+1),u_series(:,k+1),t_series(:,k+1),param)
    H(x_series(:,k+4),u_series(:,k+4),t_series(:,k+4),param)*F(x_series(:,k+4),u_series(:,k+4),t_series(:,k+4),param)*F(x_series(:,k+3),u_series(:,k+3),t_series(:,k+3),param)*F(x_series(:,k+2),u_series(:,k+2),t_series(:,k+2),param)*F(x_series(:,k+1),u_series(:,k+1),t_series(:,k+1),param)]

rank(OB)