function IEKF(sigma_discr, sigma_z, sigma_u,sigma_w_x,sigma_v_x)
global t  Delta_t
global u_meas z_meas
global mu_x Sigma2_x 
global mu_x_pred Sigma2_x_pred Q R
global param
    
    %---Begin--- One step of Kalman filter (k-th)
    
    % Prediction
    u_meas=get_u();
    % Effect of noise in u on model covariance (f_u_ is dependent on the state)
    f_u_=f_u(mu_x,u_meas,t,param);
    Q_u= f_u_*diag(sigma_u.^2)*f_u_';
    
    % Effect of noise source w_x_ on model covariance (f_w_x_ is dependent on the state)
    f_w_x_=f_w_x(mu_x,u_meas,t,param);
    Q_w_x= f_w_x_*diag(sigma_w_x.^2)*f_w_x_';
    
    % Discretization error effect on model covariance
    Q_discr=diag((sigma_discr).^2);
    
    % Error sources are assumed independent
    Q = Q_discr+Q_u+Q_w_x;
    
    f_x_=f_x(mu_x,u_meas,t,param);
    
    mu_x_pred=f(mu_x,u_meas,t,param);
    Sigma2_x_pred=f_x_*Sigma2_x*f_x_'+Q;
    I_pred=inv(Sigma2_x_pred);
    i_pred=I_pred*mu_x_pred;
    
    % Read Input and Sensor (k+1)
    t_prev=t;
    t=t+Delta_t;
    
    z_meas=get_z();
    
    % Projection / Observation
    h_u_= h_u(mu_x_pred,u_meas,t,param);
    % effect of noise in u on sensor equation covariance
    R_u= h_u_*diag(sigma_u.^2)*h_u_';
    
    h_v_x_= h_v_x(mu_x_pred,u_meas,t,param) ;
    % effect of noise source w_x_ on sensor equation covariance
    R_v_x= h_v_x_*diag(sigma_v_x.^2)*h_v_x_';
    
    R_z=diag(sigma_z.^2);
    % Sensor error, input error, and other noise sources (w_x_,...) error are assumed independent
    R=R_z+R_u+R_v_x;
    
    h_x_=h_x(mu_x_pred,u_meas,t,param);
    i_proj=h_x_'*inv(R)*(z_meas-(h(mu_x_pred,u_meas,t,param)-h_x_*mu_x_pred));
    I_proj=h_x_'*inv(R)*h_x_;
    
    % Fussion / Information Weighted Average
    I=I_pred+I_proj;
    %mu_x=inv(I)*(I_pred*mu_x_pred+I_proj*mu_x_proj);
    mu_x=inv(I)*(i_pred+i_proj);
    Sigma2_x=inv(I);
    
    %%---End--- One step of Kalman filter (k-th)

end