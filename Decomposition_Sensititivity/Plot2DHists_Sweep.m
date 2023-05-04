clc;clear;
clim = [-5 -1.5];%[0 0.05];
entropy_calc = false;

warning off;
% sweep = '22-1215-Wavy_Sweep';
% sweep = '23-0119-Wavy_Prd_Sweep';
% sweep = '23-0202-Wavy-Thick_Amp_Sweep';
% sweep = '23-0202-Wavy-Thick_Prd_Sweep';
% sweep = '23-0217-Wavy-Thick_Width_Sweep';
sweep = '23-0405-SingleCompression';

load([sweep '\sensitivity.mat'])

switch sweep 
    case {'22-1215-Wavy_Sweep','23-0202-Wavy-Thick_Amp_Sweep'}
        sin_sweep = linspace(0,2,11);
        lbl = sin_sweep;
    case {'23-0119-Wavy_Prd_Sweep','23-0202-Wavy-Thick_Prd_Sweep'}
        sin_sweep = 1:10;
        lbl = sin_sweep;
    case {'23-0217-Wavy-Thick_Width_Sweep'}
        width_sweep = 5:5:60;
        lbl = width_sweep;
    case {'23-0405-SingleCompression'}
        lbl = 0;
end

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
    colormap turbo;

    switch sweep 
        case {'22-1215-Wavy_Sweep','23-0202-Wavy-Thick_Amp_Sweep'}
            title(['Shear Wavy Simulation (' num2str(lbl(i)) ' mm Amp; 6.25mm Disp)'],'interpreter','latex')
        case {'23-0119-Wavy_Prd_Sweep','23-0202-Wavy-Thick_Prd_Sweep'}
            title(['Shear Wavy Simulation (' num2str(lbl(i)) ' Periods; 0.6mm Amp; 6.25mm Disp)'],'interpreter','latex')
        case {'23-0217-Wavy-Thick_Width_Sweep'}
            title(['Shear Wavy Simulation (' num2str(lbl(i)) 'mm Width; 4 Periods; 0.6mm Amp; 6.25mm Disp)'],'interpreter','latex')
        case {'23-0405-SingleCompression'}
            title('Biaxial Simulation - 0.15mm Disp','interpreter','latex')
    end
    if ~exist([sweep '\Pics'], 'dir')
        mkdir([sweep '\Pics'])
    end
    
    saveas(gcf,[sweep '\Pics\' sweep '_' num2str(lbl(i)) '.png']);

    if entropy_calc
        % Binwise entropy
        p = h'./sum(h(:));
        H_pwise = -sum(p.*log(p), 'all','omitnan');
        all_pwise_ent(i) = H_pwise;

        % Entropy of Gaussian
        GMModel = fitgmdist([lam{i}' k{i}'],1);
        D = 2; % 2 independent variables --> dim(X) = 2
        Sigma(:,:,i) = GMModel.Sigma;
        Mu = GMModel.mu;
        diff_entropy = (1/2)*(log(det(Sigma(:,:,i))) + (D/2)*(1 + log(2*pi)));
        all_diff_ent(i) = diff_entropy;
    end
end

close all

if entropy_calc
    figure()
    switch sweep
        case {'22-1215-Wavy_Sweep','23-0202-Wavy-Thick_Amp_Sweep'}
            ttl = 'Amplitude Sweep Entropy Plots';
            xlbl = 'Amplitude [mm]';
        case {'23-0119-Wavy_Prd_Sweep','23-0202-Wavy-Thick_Prd_Sweep'}
            ttl = 'Period Sweep Entropy Plots';
            xlbl = 'Periods [-]';
        case {'23-0217-Wavy-Thick_Width_Sweep'}
            ttl = 'Width Sweep Entropy Plots';
            xlbl = 'Width [mm]';
    end
    subplot(1,2,1)
    plot(lbl,squeeze(Sigma(2,2,:))','Color',[94, 201, 98]/255,'linewidth',2);
    sgtitle(ttl,'interpreter','latex');
    xlabel(xlbl,'interpreter','latex');
    ax = gca;
    ax.YColor = [94, 201, 98]/255;
    addaxis(lbl,squeeze(Sigma(1,1,:))','Color',[33, 145, 140]/255,'linewidth',2)
    legend({'Gaussian Std. Dev., $\sigma_k$','Gaussian Std. Dev., $\sigma_\lambda$'},...
        'interpreter','latex','location','best')

    subplot(1,2,2)
    plot(lbl,all_diff_ent,'Color',[59, 82, 139]/255,'linewidth',2)
    xlabel(xlbl,'interpreter','latex');
    ax = gca;
    ax.YColor = [59, 82, 139]/255;
    addaxis(lbl,all_pwise_ent,'Color',[68, 1, 84]/255,'linewidth',2)
    legend({'Entropy of Gaussian, $p_{gauss}$','Binwise Entropy, $p_{bin}$'},...
        'interpreter','latex','location','best')
%     addaxislabel(1,'Gaussian Std. Dev., \sigma_k','interpreter','latex');
%     addaxislabel(2,'Entropy of Gaussian, p_{gauss}','interpreter','latex');
%     addaxislabel(3,'Binwise Entropy, p_{bin}','interpreter','latex');

    saveas(gcf,[sweep '\Pics\' sweep '_Gaussian.png']);
end