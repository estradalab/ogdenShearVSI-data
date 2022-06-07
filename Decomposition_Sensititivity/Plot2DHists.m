%% Run sensitivity for S_alpha functions
sensitivity_metric_scrap
close all

type = 'Sim'; % Sim - Simulation; Exp - Experiment

%% k vs. lam Holes
switch type
    case 'Sim'
        load('22-0317-Updated-Dataset\Simulations\Holes\sensitivity.mat')
    case 'Exp'
        load('22-0317-Updated-Dataset\Experimental\Holes\sensitivity.mat')
end
m = 2;
figure(m*10); sgtitle(['Holes for ' type],'interpreter','latex')
subplot(1,3,1);
[lam_bin,k_bin,h] = ndhist(lam{1}(:),k{1}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{1,1} = lam_bin;
K_BIN{1,1} = k_bin;
H{1} = h';
% imagesc(h')
% set(gca,'YDir','normal'); 
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
% set(gca,'YDir','normal'); 
axis square;
ylabel('k')
xlabel('$\lambda$','interpreter','latex')
title('u=5$\,$mm','interpreter','latex')
colormap turbo;

saveas(gcf,['k_lam_SensitivityPlots/' type '/histo_holes.png']);
saveas(gcf,['k_lam_SensitivityPlots/' type '/histo_holes.pdf']);

%% k vs. lam No Holes
switch type
    case 'Sim'
        load('22-0317-Updated-Dataset\Simulations\No_Holes\sensitivity.mat')
    case 'Exp'
        load('22-0317-Updated-Dataset\Experimental\No_Holes\sensitivity.mat')
end
figure(m*10+1);
sgtitle(['No Holes for ' type],'interpreter','latex')
subplot(1,3,1);
[lam_bin,k_bin,h] = ndhist(lam{1}(:),k{1}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{2,1} = lam_bin;
K_BIN{2,1} = k_bin;
H{3} = h';
% imagesc(h')
% set(gca,'YDir','normal'); 
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
% set(gca,'YDir','normal'); 
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
% imagesc(h')
% set(gca,'YDir','normal'); 
axis square;
ylabel('k')
xlabel('$\lambda$','interpreter','latex')
title('u=7$\,$mm','interpreter','latex')
catch
end
colormap turbo;

saveas(gcf,['k_lam_SensitivityPlots/' type '/histo_no_holes.png']);
saveas(gcf,['k_lam_SensitivityPlots/' type '/histo_no_holes.pdf']);

%close all;

%% k vs. lam Sensitivity Plots

figure; % imagesc(squeeze(Salpha_all(:,1,:))') % a = -5
[LAM,K] = ndgrid(lam_bin,k_bin);
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,1,:)),500);
set(hand,'LineColor','none');
axis square;
title([type ': $\alpha=-5$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
saveas(gcf,['k_lam_SensitivityPlots/' type '/sens_alpha-5.png']);
saveas(gcf,['k_lam_SensitivityPlots/' type '/sens_alpha-5.pdf']);

figure; % imagesc(squeeze(Salpha_all(:,1001,:))') % a = 5
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,1001,:)),500);
set(hand,'LineColor','none');
axis square;
title([type ': $\alpha=5$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
saveas(gcf,['k_lam_SensitivityPlots/' type '/sens_alpha5.png']);
saveas(gcf,['k_lam_SensitivityPlots/' type '/sens_alpha5.pdf']);

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
        mat_temp = squeeze(Salpha_all(:,i,:));
        mat_temp = mat_temp./sum(mat_temp(:)); % Normalized
        S_lookup{test}(i,1) = alpha(i);
        S_lookup{test}(i,2) = sum(H{test}(:).*mat_temp(:));
    end
    if test == 5
        S_hat_lookup = S_lookup;
        for j = 1:5
            nrmlz(j) = sum(H{j}(:));
            S_hat_lookup{j}(:,2) = S_lookup{j}(:,2)/sum(H{j}(:)); % Normalized
            plot(S_hat_lookup{j}(:,1),S_hat_lookup{j}(:,2))
            hold on
        end
    end
end
title(['$S_{metric}$(test,$\alpha$) for ' type],'interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
ylabel('$S_{metric}=\sum_{i,j}\hat{hist}\left(k_i,\lambda_j,test\right)\;.*\;\hat{S}_\alpha\left(k_i,\lambda_j,\alpha\right)$',...
    'interpreter','latex')
legend('Holes; $U_d = 2.5 mm$','Holes; $U_d = 5 mm$','No Holes; $U_d = 2.5 mm$',...
    'No Holes; $U_d = 5 mm$', 'No Holes; $U_d = 7 mm$','interpreter','latex')
saveas(gcf,['k_lam_SensitivityPlots/' type '/sens_metric.png'])
saveas(gcf,['k_lam_SensitivityPlots/' type '/sens_metric.pdf'])

save(['S_metrics_histos_' type '.mat'],H,S_hat_lookup,S_lookup);
    
    


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
%     saveas(gcf,['k_lam_SensitivityPlots/sens_' num2str(i) '.png']);
% end

% figure; imagesc(squeeze(Salpha_all(:,501,:))') % a = 2
% set(gca,'YDir','normal'); axis image;
% axis square;
% 
% 
% figure; imagesc(squeeze(Salpha_all(:,1,:))') % a = -3
% set(gca,'YDir','normal'); axis image;
% axis square;