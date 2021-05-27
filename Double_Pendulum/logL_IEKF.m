%%
function logL=logL_IEKF(diag_Sigma_discr,diag_Sigma_z,diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x,mu_x_0,diag_Sigma_x_0)
global t t_0 t_end Delta_t
global u_meas z_meas 
global mu_x Sigma2_x
global mu_x_pred Sigma2_x_pred R
global seed
global param

%% Kalman Filter Loop

% Starting simulation time and length
logL=0;
%% Read Input and Sensor at first (k=1)

rng(seed);

t=t_0;
mu_x=mu_x_0;
Sigma2_x=diag((diag_Sigma_x_0).^2);
rng(seed); u_meas=get_u(); z_meas=get_z();

for k=1:t_end/Delta_t

    IEKF(diag_Sigma_discr, diag_Sigma_z, diag_Sigma_u,diag_Sigma_w_x,diag_Sigma_v_x);
    
    %% Log-Likelihood computation
    h_x_=h_x(mu_x_pred,u_meas,t,param);
    Sigma2_Innovation=h_x_*Sigma2_x_pred*h_x_'+R;
    Innovation=z_meas-h(mu_x_pred,u_meas,t,param);
    logL=logL+1/2*(log(det(Sigma2_Innovation))+Innovation'*inv(Sigma2_Innovation)*Innovation);
end
end