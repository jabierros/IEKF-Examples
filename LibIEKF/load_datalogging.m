function load_datalogging(sol, datalogging_string)
% Datalogged variables
global t u_meas u_actual z_meas z_actual x_actual mu_x Sigma2_x sigma_x 

load(sol,'-ascii')

column_end=0;
for i=1:size(datalogging_string,1)
    n_columns=size(eval(datalogging_string{i}),1);
    column_start=column_end+1;
    column_end=column_end+n_columns;
    assignin('caller',[datalogging_string{i},'_series'], sol(:,column_start:column_end));
end
clear sol
end