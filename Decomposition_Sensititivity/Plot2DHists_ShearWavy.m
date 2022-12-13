%% Run sensitivity for S_alpha functions
type = {'Sim','Exp'};
flt = 'Filter'; % 'Filter','Filter_None'
fn = 'cax-5to0LOG';
clim = [-5 -1.5];%[0 0.05];

warning off;

load('22-1212-Shear_Wavy\sensitivity.mat')

m = 2;
figure(m*10); sgtitle('Shear Wavy Simulation (6.5 mm Disp)','interpreter','latex')
subplot(2,2,1);
[lam_bin,k_bin,h] = ndhist(lam{1}(:),k{1}(:),'axis',[1,0;2,1]-0.005);
LAM_BIN{1,1} = lam_bin;
K_BIN{1,1} = k_bin;
H{1} = h';
imagesc([1 2],[0 1],h/sum(h(:)));
% caxis(clim); 
colormap turbo;
set(gca,'YDir','normal'); 

axis square;
ylabel('k','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')
title('Rectangular','interpreter','latex')
colormap turbo; 


subplot(2,2,2);
[lam_bin,k_bin,h] = ndhist(lam{2}(:),k{2}(:),'axis',[1,0;2,1]-0.005);
LAM_BIN{1,2} = lam_bin;
K_BIN{1,2} = k_bin;
H{2} = h';
imagesc([1 2],[0 1],h/sum(h(:))); c = colorbar;
c.Label.String = 'h/sum(h)';
% caxis(clim);
colormap turbo;
set(gca,'YDir','normal'); 

axis square;
ylabel('k','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')
title('Wavy','interpreter','latex')
colormap turbo;

subplot(2,2,1); caxis(c.Limits); c = colorbar;
c.Label.String = 'h/sum(h)';

subplot(2,2,3);
[lam_bin,k_bin,h] = ndhist(lam{1}(:),k{1}(:),'axis',[1,0;2,1]-0.005);
LAM_BIN{1,1} = lam_bin;
K_BIN{1,1} = k_bin;
H{1} = h';
imagesc([1 2],[0 1],log10(h/sum(h(:))));
caxis(clim); c = colorbar; c.Label.String = 'log(h/sum(h))';
colormap turbo;
set(gca,'YDir','normal'); 

axis square;
ylabel('k','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')
title('Rectangular - log','interpreter','latex')
colormap turbo; 


subplot(2,2,4);
[lam_bin,k_bin,h] = ndhist(lam{2}(:),k{2}(:),'axis',[1,0;2,1]-0.005);
LAM_BIN{1,2} = lam_bin;
K_BIN{1,2} = k_bin;
H{2} = h';
imagesc([1 2],[0 1],log10(h/sum(h(:))))
caxis(clim); c = colorbar; c.Label.String = 'log(h/sum(h))';
colormap turbo;
set(gca,'YDir','normal'); 

axis square;
ylabel('k','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')
title('Wavy - log','interpreter','latex')
colormap turbo;

saveas(gcf,'22-1212-Shear_Wavy\2D_Hist.png');