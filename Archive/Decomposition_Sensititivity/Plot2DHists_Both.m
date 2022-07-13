%% Run sensitivity for S_alpha functions
sensitivity_metric_scrap
close all

type = {'Sim','Exp'};
flt = 'Filter'; % 'Filter','Filter_None'
fn = 'cax-5to0LOG';
clim = [-5 0];%[0 0.05];

for ii = 1:2
%% k vs. lam Holes
switch flt
    case 'Filter'
        flt_str = '';
    case 'Filter_None'
        flt_str = '_No_Filter';
end

switch type{ii}
    case 'Sim'
        load(['22-0317-Updated-Dataset\Simulations\Holes\sensitivity.mat'])
    case 'Exp'
        load(['22-0317-Updated-Dataset\Experimental\Holes' flt_str '\sensitivity.mat'])
end
m = 2;
figure(m*10); sgtitle(['Holes for ' type{ii}],'interpreter','latex')
subplot(1,3,1);
[lam_bin,k_bin,h] = ndhist(lam{1}(:),k{1}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{1,1} = lam_bin;
K_BIN{1,1} = k_bin;
H{1} = h';
% imagesc(h')
%imagesc(h/sum(h(:))); colormap; caxis(clim); 
imagesc(log10(h/sum(h(:)))); caxis(clim); 
colormap turbo;
set(gca,'YDir','normal'); 

axis square;
ylabel('k')
xlabel('$\lambda$','interpreter','latex')
title('u=2.5$\,$mm','interpreter','latex')
colormap turbo; 

% subplot(1,2,2);
% contourf(h',15,'LineStyle','none');
% set(gca,'YDir','normal');


subplot(1,3,2);
[lam_bin,k_bin,h] = ndhist(lam{2}(:),k{2}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{1,2} = lam_bin;
K_BIN{1,2} = k_bin;
H{2} = h';
% imagesc(h')
%imagesc(h/sum(h(:))); colormap; caxis(clim); 
imagesc(log10(h/sum(h(:)))); caxis(clim); colormap turbo;
set(gca,'YDir','normal'); 

axis square;
ylabel('k')
xlabel('$\lambda$','interpreter','latex')
title('u=5$\,$mm','interpreter','latex')
colormap turbo;

saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/histo_holes' fn '.png']);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/histo_holes' fn '.pdf']);

%% k vs. lam No Holes
switch type{ii}
    case 'Sim'
        load(['22-0317-Updated-Dataset\Simulations\No_Holes\sensitivity.mat'])
    case 'Exp'
        load(['22-0317-Updated-Dataset\Experimental\No_Holes' flt_str '\sensitivity.mat'])
end
figure(m*10+1);
sgtitle(['No Holes for ' type{ii}],'interpreter','latex')
subplot(1,3,1);
[lam_bin,k_bin,h] = ndhist(lam{1}(:),k{1}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{2,1} = lam_bin;
K_BIN{2,1} = k_bin;
H{3} = h';
%imagesc(h/sum(h(:))); colormap; caxis(clim); 
imagesc(log10(h/sum(h(:)))); caxis(clim);
colormap turbo;
% imagesc(h')
set(gca,'YDir','normal'); 

axis square;
ylabel('k')
xlabel('$\lambda$','interpreter','latex')
title('u=2.5$\,$mm','interpreter','latex')
colormap turbo; 

subplot(1,3,2);
[lam_bin,k_bin,h] = ndhist(lam{2}(:),k{2}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{2,2} = lam_bin;
K_BIN{2,2} = k_bin;
H{4} = h';
% imagesc(h')
%imagesc(h/sum(h(:))); colormap; caxis(clim); 
imagesc(log10(h/sum(h(:)))); caxis(clim);
colormap turbo;
set(gca,'YDir','normal');

axis square;
ylabel('k')
xlabel('$\lambda$','interpreter','latex')
title('u=5$\,$mm','interpreter','latex')
colormap turbo;

try
subplot(1,3,3);
[lam_bin,k_bin,h] = ndhist(lam{3}(:),k{3}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{2,3} = lam_bin;
K_BIN{2,3} = k_bin;
H{5} = h';
%imagesc(h/sum(h(:))); colormap; caxis(clim); 
imagesc(log10(h/sum(h(:)))); caxis(clim);
colormap turbo;
set(gca,'YDir','normal');
% imagesc(h')
% set(gca,'YDir','normal'); 
axis square;
ylabel('k')
xlabel('$\lambda$','interpreter','latex')
title('u=7$\,$mm','interpreter','latex')
catch
end
colormap turbo;

saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/histo_no_holes' fn '.png']);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/histo_no_holes' fn '.pdf']);

%close all;

%% k vs. lam Sensitivity Plots

figure; % imagesc(squeeze(Salpha_all(:,1,:))') % a = -5
[LAM,K] = ndgrid(lam_bin,k_bin);
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,1,:)),50);
set(hand,'LineColor','none');
axis square;
title([type{ii} ': $\alpha=-5$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
colormap magma;
set(gcf,'renderer','painters');
caxis([0 8]);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha-5.png']);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha-5.pdf']);

figure; % imagesc(squeeze(Salpha_all(:,1,:))') % a = -3
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,201,:)),50);
set(hand,'LineColor','none');
axis square;
title([type{ii} ': $\alpha=-3$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
colormap magma;
set(gcf,'renderer','painters');
caxis([0 8]);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha-3.png']);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha-3.pdf']);

figure; % imagesc(squeeze(Salpha_all(:,1,:))') % a = -1
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,401,:)),50);
set(hand,'LineColor','none');
axis square;
title([type{ii} ': $\alpha=-1$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
colormap magma;
set(gcf,'renderer','painters');
caxis([0 8]);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha-1.png']);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha-1.pdf']);

figure; % imagesc(squeeze(Salpha_all(:,1,:))') % a = 1
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,601,:)),50);
set(hand,'LineColor','none');
axis square;
title([type{ii} ': $\alpha=1$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
colormap magma;
set(gcf,'renderer','painters');
caxis([0 8]);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha1.png']);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha1.pdf']);

figure; % imagesc(squeeze(Salpha_all(:,1,:))') % a = 3
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,801,:)),50);
set(hand,'LineColor','none');
axis square;
title([type{ii} ': $\alpha=3$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
colormap magma;
set(gcf,'renderer','painters');
caxis([0 8]);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha3.png']);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha3.pdf']);

figure; % imagesc(squeeze(Salpha_all(:,1001,:))') % a = 5
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,1001,:)),50);
set(hand,'LineColor','none');
axis square;
title([type{ii} ': $\alpha=5$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
colormap magma;
set(gcf,'renderer','painters');
caxis([0 8]);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha5.png']);
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_alpha5.pdf']);

%% S_alpha{test}(metric;alpha) Plots
% Tests: 
% 1 - Holes; 2.5 mm
% 2 - Holes; 5 mm
% 3 - No Holes; 2.5 mm
% 4 - No Holes; 5 mm
% 5 - No Holes; 7 mm
figure();
for test = 1:5
    for i = 1:length(alpha)
        mat_temp{i} = squeeze(Salpha_all(:,i,:));
        mat_temp2{i} = squeeze(Smu_all(:,i,:));
        % mat_temp{i} = mat_temp{i}./sum(mat_temp{i}(:)); % Normalized
        S_lookup{test}(i,1) = alpha(i);
        S_lookup{test}(i,2) = sum(H{test}(:).*mat_temp{i}(:));
        S_lookup{test}(i,3) = sum(H{test}(:).*mat_temp2{i}(:));
    end
    if test == 5
        S_hat_lookup = S_lookup;
        for j = 1:5
            nrmlz(j) = sum(H{j}(:));
            S_hat_lookup{j}(:,2) = S_lookup{j}(:,2)/sum(H{j}(:)); % Normalized
            S_hat_lookup{j}(:,3) = S_lookup{j}(:,3)/sum(H{j}(:));
            plot(S_hat_lookup{j}(:,1),S_hat_lookup{j}(:,2))
            hold on
        end
    end
end
title(['$S_{metric}$(test,$\alpha$) for ' type{ii}],'interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
ylabel('$S_{metric}=\sum_{i,j}\hat{hist}\left(k_i,\lambda_j,test\right)\;.*\;\hat{S}_\alpha\left(k_i,\lambda_j,\alpha\right)$',...
    'interpreter','latex')
legend('Holes; $U_d = 2.5 mm$','Holes; $U_d = 5 mm$','No Holes; $U_d = 2.5 mm$',...
    'No Holes; $U_d = 5 mm$', 'No Holes; $U_d = 7 mm$','interpreter','latex')
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_metric.png'])
saveas(gcf,[flt '/k_lam_SensitivityPlots/' type{ii} '/sens_metric.pdf'])

Salpha_ALL = mat_temp;
Smu_ALL = mat_temp2;
save([flt '/k_lam_SensitivityPlots/S_metrics_histos_' type{ii} '.mat'],'H','S_hat_lookup','S_lookup','Salpha_ALL');

H_{ii} = H; S_hat_{ii} = S_hat_lookup; S_{ii} = S_lookup; Salpha_all_{ii} = Salpha_ALL; Smu_all_{ii} = Smu_ALL;
    
end 

%% S_alpha{test}(metric;alpha) Plots (Both sim and exp)
save([flt '/k_lam_SensitivityPlots/S_metrics_histos_both.mat'],'H_','S_hat_','S_','Salpha_all_','Smu_ALL');
figure()
for ii = 1:length(H_)
    if ii == 1
        for jj = 1:length(H_{ii})
            plot(S_hat_{ii}{jj}(:,1),S_hat_{ii}{jj}(:,2))
            hold on
        end
    else
        for jj = 1:length(H_{ii})
            plot(S_hat_{ii}{jj}(:,1),S_hat_{ii}{jj}(:,2))
            if jj == 1
                set(gca,'ColorOrderIndex',1)
            end
        end
    end
end
title(['$S_{metric}$(test,$\alpha$)'],'interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
ylabel('$S_{metric}=\sum_{i,j}\hat{hist}\left(k_i,\lambda_j,test\right)\;.*\;\hat{S}_\alpha\left(k_i,\lambda_j,\alpha\right)$',...
    'interpreter','latex')
legend('Sim; Holes; $U_d = 2.5 mm$','Sim; Holes; $U_d = 5 mm$','Sim; No Holes; $U_d = 2.5 mm$',...
    'Sim; No Holes; $U_d = 5 mm$', 'Sim; No Holes; $U_d = 7 mm$','Exp; Holes; $U_d = 2.5 mm$',...
    'Exp; Holes; $U_d = 5 mm$','Exp; No Holes; $U_d = 2.5 mm$',...
    'Exp; No Holes; $U_d = 5 mm$', 'Exp; No Holes; $U_d = 7 mm$','interpreter','latex')

saveas(gcf,[flt '/k_lam_SensitivityPlots/sens_metric_both.png'])
% saveas(gcf,[flt '/k_lam_SensitivityPlots/sens_metric_both.pdf'])

%% Other

% a = [-4.94,-4.22,-3.87,-2,-1.84,-1.33,1,1.45,1.5,1.65,1.89,1.91,2,2.86,3.03,3.3];
% 
% for i = 1:length(a)
%     idx = round((a(i) + 5.01)*100);
%     figure; imagesc(squeeze(Salpha_all(:,idx,:))') % a = 2
%     set(gca,'YDir','normal'); axis image;
%     title(['$\alpha=$',num2str(a(i))],'interpreter','latex')
%     ylabel('k'); xlabel('$\lambda$','interpreter','latex')
%     axis square;
%     saveas(gcf,[flt '/k_lam_SensitivityPlots/sens_' num2str(i) '.png']);
% end

% figure; imagesc(squeeze(Salpha_all(:,501,:))') % a = 2
% set(gca,'YDir','normal'); axis image;
% axis square;
% 
% 
% figure; imagesc(squeeze(Salpha_all(:,1,:))') % a = -3
% set(gca,'YDir','normal'); axis image;
% axis square;