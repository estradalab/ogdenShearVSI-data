clc;clear;
%% Run sensitivity for S_alpha functions
type = {'Sim','Exp'};
flt = 'Filter'; % 'Filter','Filter_None'
fn = 'cax-5to0LOG';
clim = [-5 -1.5];%[0 0.05];

warning off;

load('22-1215-Wavy_Sweep\sensitivity.mat')
sin_sweep = linspace(0,2,11);

for i = 1:length(k)
    figure(i);
    [lam_bin,k_bin,h] = ndhist(lam{i}(:),k{i}(:),'axis',[1,0;2,1]-0.005);
    LAM_BIN{i} = lam_bin;
    K_BIN{i} = k_bin;
    H{i} = h';
    imagesc([1 2],[0 1],log10(h/sum(h(:))));
    caxis(clim); c = colorbar; c.Label.String = 'log(h/sum(h))';
    colormap turbo;
    set(gca,'YDir','normal'); 
    
    axis square;
    ylabel('k','interpreter','latex')
    xlabel('$\lambda$','interpreter','latex')
     title(['Shear Wavy Simulation (' num2str(sin_sweep(i)) ' mm Amp; 6.25mm Disp)'],'interpreter','latex')
    colormap turbo;

    saveas(gcf,['22-1215-Wavy_Sweep\Pics\2D_Hist_Amp_' num2str(sin_sweep(i)) '.png']);
end