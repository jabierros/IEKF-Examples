function Plotting(fig_dir, datalogging_string)
global t_0 t_end
global x_actual_0
global param
%global n_x n_u n_z
global x_string u_string z_string

[mkdir_success,mkdir_message]=mkdir(fig_dir);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

fontsize=12;

load_datalogging('sol.dat', datalogging_string);

mu_x_error_series=mu_x_series-x_actual_series;
%--------------------------------------------------------------------------
figIdx=1;

fig=figure(figIdx);
figname='x_actual';
plot(t_series,x_actual_series);
set(gca, 'YScale', 'linear');
legend_list=[];
n_x=size(x_actual_series,2);
for i=1:n_x
    legend_list{i}=strcat('$',x_string(i),'^{tr}$');
end
legH=legend(legend_list);
titH=title('$\mathbf{x}^{tr}$ (Actual State)');
set(fig,'units','normalized');
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize)
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);
%--------------------------------------------------------------------------
x_discr_series=[];
x_=x_actual_0;
for k=1:length(t_series)
 x_discr_series=[x_discr_series,x_];
 x_=f(x_, u_actual_series(k,:)', t_series(k),param);
end

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='x_discr';
plot(t_series,x_discr_series);
set(gca, 'YScale', 'linear');
legend_list=[];
for i=1:n_x
    legend_list{i}=strcat('$',x_string(i),'^{discr}$');
end
legH=legend(legend_list);
titH=title('$\mathbf{x}^{discr}$ (integrated with discr.)');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);
%--------------------------------------------------------------------------
pred_error_series=zeros(size(x_actual_series,2),size(t_series,1));
for k=1:size(t_series,1)-1
    pred_error_series(:,k+1)=abs((f(x_actual_series(k,:)',u_actual_series(k,:)',t_series(k,:)',param)-x_actual_series(k,:)')-(x_actual_series(k+1,:)'-x_actual_series(k,:)'));
end
pred_error_series=pred_error_series';

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='discr_error';
h=plot(t_series,pred_error_series);
set(gca, 'YScale', 'log');
legend_list=[];
for i=1:n_x
    legend_list{i}=strcat('$\varepsilon_{',x_string(i),'}$');
end
legH=legend(legend_list);
titH=title('Discretization error ($\left|\mathbf{f}(\mathbf{x}^{tr}_k)-\mathbf{x}^{tr}_{k+1}\right|$)');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
%--------------------------------------------------------------------------
dcm_obj = datacursormode(gcf);
for i=1:length(h)
[y_max,x_max]=max(pred_error_series(:,i));
% hDatatip(i)=createDatatip(dcm_obj, h(i), [t_series(x_max),pred_error_series(x_max,i)]);
% set(hDatatip(i),'Position',[t_series(x_max),pred_error_series(x_max,i)]);
% set(hDatatip(i),'interpreter','latex');
end

set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='u_meas_u_actual';
xlim([t_0,t_end])
n_u=size(u_meas_series,2);
if n_u>0
plot(t_series,u_meas_series,'-',t_series,u_actual_series,'--');
end
set(gca, 'YScale', 'linear');
legend_list=[];
for i=1:n_u
    legend_list{i}=strcat('$\mathrm{meas}(',u_string(i),')$');
end
for i=1:n_u
    legend_list{i+n_u}=strcat('$',u_string(i),'$');
end
legH=legend(legend_list);
titH=title('$\mathrm{meas}(\mathbf{u})$ vs $\mathbf{u}$ ');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);
%--------------------------------------------------------------------------
figIdx = figIdx + 1;
fig=figure(figIdx);
figname='u_meas_minus_u_actual';
legend_list=[];
if n_u>0
plot(t_series,u_meas_series-u_actual_series,'-');
xlim([t_0,t_end])
end
set(gca, 'YScale', 'linear');
legend_list=[];
for i=1:n_u
    legend_list{i}=strcat('$\mathrm{meas}(',u_string(i),')-',u_string(i),'$');
end
legH=legend(legend_list);
titH=title('$\mathrm{meas}(\mathbf{u})-\mathbf{u}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);
%--------------------------------------------------------------------------
figIdx = figIdx + 1;
fig=figure(figIdx);
figname='z_meas_z_actual';
plot(t_series,z_meas_series,'-',t_series,z_actual_series,'--');
set(gca, 'YScale', 'linear');
legend_list=[];
n_z=size(z_meas_series,2);
for i=1:n_z
    legend_list{i}=strcat('$\mathrm{meas}(',z_string(i),')$');
end
for i=1:n_z
    legend_list{i+n_z}=strcat('$',z_string(i),'$');
end
legH=legend(legend_list);
titH=title('$\mathrm{meas}(\mathbf{z})$ vs $\mathbf{z}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);
%--------------------------------------------------------------------------
figIdx = figIdx + 1;
fig=figure(figIdx);
figname='z_meas_minus_z_actual';
plot(t_series,z_meas_series-z_actual_series,'-');
set(gca, 'YScale', 'linear');
legend_list=[];
for i=1:n_z
    legend_list{i}=strcat('$\mathrm{meas}(',z_string(i),')-',z_string(i),'$');
end
legH=legend(legend_list);
titH=title('$\mathrm{meas}(\mathbf{z})-\mathbf{z}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);
%--------------------------------------------------------------------------
figIdx = figIdx + 1;
fig=figure(figIdx);
figname='mu_x_x_actual';
plot(t_series,mu_x_series,'-',t_series,x_actual_series,'--');
set(gca, 'YScale', 'linear');
legend_list=[];
for i=1:n_x
    legend_list{i}=strcat('$\hat{\mu}_{',x_string(i),'}$');
end
for i=1:n_x
    legend_list{i+n_x}=strcat('$',x_string(i),'$');
end
legH=legend(legend_list);
xlim([0 t_end]);
% ylim([-6 6]);
titH=title('$\mathbf{x}$ vs $\hat{{\mu}}_{\mathbf{x}}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);
%--------------------------------------------------------------------------
figIdx = figIdx + 1;
fig=figure(figIdx);
set(fig,'units','normalized'); 
plot(t_series,mu_x_series-x_actual_series,'-');
set(gca, 'YScale', 'linear');
legend_list=[];
for i=1:n_x
    legend_list{i}=strcat('$\hat{\mu}_{',x_string(i),'}-',x_string(i),'$');
end
legH=legend(legend_list);
xlim([0 t_end]);
% ylim([-1 1]);
set(legH,'interpreter','latex','Fontsize',fontsize);
titH=title('$\hat{\mu}_{\mathbf{x}}-\mathbf{x}$');
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
figname='mu_x_minus_x_actual';
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);
%--------------------------------------------------------------------------
% Real filter error statistics @ lim k -> infty
num_samples_statistic=100;

lim_mu_x_error=mu_x_error_series(end-num_samples_statistic:end,:);
lim_mu_x_error_mean=mean(lim_mu_x_error);
lim_mu_x_error_std=std(lim_mu_x_error);
sqrt_lim_mu_x_error_squared_mean=mean((lim_mu_x_error).^2).^0.5;

figIdx = figIdx + 1;
fig=figure(figIdx);
figname='sigma_x';
plot(t_series,sigma_x_series,'-',t_series,ones(size(sigma_x_series,1),1)*sqrt_lim_mu_x_error_squared_mean,'--');
set(gca, 'YScale', 'log');
legend_list=[];
for i=1:n_x
    legend_list{i}=strcat('${\sigma}_{',x_string(i),'}$');
end
for i=1:n_x
    legend_list{i+n_x}=strcat('${\mathrm{mean}((\hat{\mu}_{',x_string(i),'}-{',x_string(i),'})^2)^{\frac{1}{2}}}$');
end
legH=legend(legend_list);
yl = ylim;
%ylim([1e-3,yl(2)]);
titH=title('$\mathrm{diag}({{\Sigma}}^2_{\mathbf{x}})^{\frac{1}{2}}$');
set(fig,'units','normalized'); 
set(legH,'interpreter','latex','Fontsize',fontsize);
set(titH,'interpreter','latex','Fontsize',fontsize);
set(fig.CurrentAxes,'FontSize',fontsize);
set(gcf,'renderer','painters');saveas(gcf,[fig_dir,'/',figname],'png');
%system( ['cd ',fig_dir,' ; epstopdf ',figname,'.eps ; cd ..']);
%--------------------------------------------------------------------------
end