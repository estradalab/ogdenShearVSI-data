clc;clear;
selpath = uigetdir;
[~,folderName] = fileparts(selpath);
load([selpath '\output.mat']);

maxLam = 2;
bin_res = 0.01;
clim = [-5 -1.5];

decomp = questdlg('Select decompositon choice','Decomposition',...
    'K_Lambda','Criscione','Criscione');

switch decomp
    case 'K_Lambda'
        figure()
        [~,~,h] = ndhist(output.lam{1}(:),output.k{1}(:),'axis',[1,0;maxLam,1],'fixed_bin_res',bin_res);
        imagesc([1 maxLam],[0 1],log10(h/sum(h(:))));
        axis square
        colormap turbo;
        set(gca,'YDir','normal'); 
        ylabel('k, Stretch state','interpreter','latex')
        xlabel('$\lambda$, Stretch amplitude','interpreter','latex')
        title(['K vs Lambda Decomposition ' folderName],'interpreter','latex')
        c = colorbar; caxis(clim);
        c.Label.String = 'log(h/sum(h(:)))';
        saveas(gcf,[selpath '\2Dhist_k_lam.png']); saveas(gcf,[selpath '\2Dhist_k_lam.pdf'])
    case 'Criscione'
        k_2 = 0:bin_res:floor(log(maxLam)/sqrt(2/3)/bin_res)*bin_res;
        k_3=-1:2*bin_res:1;
        figure()
        [~,~,h] = ndhist(output.k2{1}(:),output.k3{1}(:),'axis',[min(k_2) max(k_2) min(k_3) max(k_3)],'fixed_bin_res',bin_res,'scale_y',2);
        imagesc([min(k_2) max(k_2)],[min(k_3) max(k_3)],log10(h/sum(h(:))));
        axis square
        colormap turbo;
        set(gca,'YDir','normal'); 
        ylabel('$k_3$, Strain state','interpreter','latex')
        xlabel('$k_2$, Log strain amplitude','interpreter','latex')
        title(['k2 vs k3 Decomposition ' folderName],'interpreter','latex')
        c = colorbar; caxis(clim);
        c.Label.String = 'log(h/sum(h(:)))';
        saveas(gcf,[selpath '\2Dhist_k2_k3.png']); saveas(gcf,[selpath '\2Dhist_k2_k3.pdf'])
end