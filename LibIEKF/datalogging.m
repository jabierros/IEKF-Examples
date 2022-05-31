function datalogging(fid, datalogging_string)
% Datalogged variables
global t u_meas u_actual z_meas z_actual x_actual mu_x Sigma2_x sigma_x

sigma_x=diag(Sigma2_x).^0.5;

for i=1:size(datalogging_string,1)
    var=eval(datalogging_string{i});
    for j=1:size(var,1)
     fprintf(fid,'%d ',var(j,1));
    end
end
fprintf(fid,'\n');
end