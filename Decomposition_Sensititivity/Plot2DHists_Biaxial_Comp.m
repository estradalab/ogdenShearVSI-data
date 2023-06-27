%% Run sensitivity for S_alpha functions
sensitivity_metric_scrap
close all

type = {'Exp','Sim'};
flt = 'Filter'; % 'Filter','Filter_None'
fn = 'cax-5to0LOG';
clim = [-5 -1.5];%[0 0.05];

warning off;
for ii = 1:1 % Change to 2 when doing it with simulation
%% k vs. lam (1)Biaxial/(2)Uniaxial/(3)MR-u Uniaxial

load('23-0626-Biaxial_Uniaxial_Comp\Experimental\sensitivity_biaxial.mat')

figure(1);
[lam_bin,k_bin,h] = ndhist(lam{1}(:),k{1}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{1} = lam_bin;
K_BIN{1} = k_bin;
H{1} = h';
imagesc([1 1.8],[0 1],log10(h/sum(h(:))));
caxis(clim); c = colorbar; c.Label.String = 'log(h/sum(h))';
colormap turbo;
set(gca,'YDir','normal'); 

axis square;
ylabel('k','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')

title('Biaxial stretch (4.5 mm + constrained)','interpreter','latex')

saveas(gcf,'23-0626-Biaxial_Uniaxial_Comp/biaxial_2dhist.png');
saveas(gcf,'23-0626-Biaxial_Uniaxial_Comp/biaxial_2dhist.pdf');

load('23-0626-Biaxial_Uniaxial_Comp\Experimental\sensitivity_uniaxial.mat')

figure(2);
[lam_bin,k_bin,h] = ndhist(lam{1}(:),k{1}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{2} = lam_bin;
K_BIN{2} = k_bin;
H{2} = h';
imagesc([1 1.8],[0 1],log10(h/sum(h(:))));
caxis(clim); c = colorbar; c.Label.String = 'log(h/sum(h))';
colormap turbo;
set(gca,'YDir','normal'); 

axis square;
ylabel('k','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')

title('Uniaxial stretch (1 mm)','interpreter','latex')

saveas(gcf,'23-0626-Biaxial_Uniaxial_Comp/uniaxial_2dhist.png');
saveas(gcf,'23-0626-Biaxial_Uniaxial_Comp/uniaxial_2dhist.pdf');

load('22-0325-Ogden_Uniaxial\Experimental\sensitivity.mat')

figure(3);
[lam_bin,k_bin,h] = ndhist(lam{10}(:),k{10}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN{3} = lam_bin;
K_BIN{3} = k_bin;
H{3} = h';
imagesc([1 1.8],[0 1],log10(h/sum(h(:))));
caxis(clim); c = colorbar; c.Label.String = 'log(h/sum(h))';
colormap turbo;
set(gca,'YDir','normal'); 

axis square;
ylabel('k','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')

title('MR-u paper uniaxial stretch (7 mm)','interpreter','latex')

saveas(gcf,'23-0626-Biaxial_Uniaxial_Comp/mr_u_2dhist.png');
saveas(gcf,'23-0626-Biaxial_Uniaxial_Comp/mr_u_2dhist.pdf');

%% S_alpha{test}(metric;alpha) Plots
% Tests: 
% 1 - Biaxial; 4.5 mm
% 2 - Uniaxial; 1 mm
% 3 - MR-u Uniaxial; 7 mm
figure();
for test = 1:3
    for i = 1:length(alpha)
        mat_temp{i} = squeeze(Salpha_all(:,i,:));
        mat_temp2{i} = squeeze(Smu_all(:,i,:));
        % mat_temp{i} = mat_temp{i}./sum(mat_temp{i}(:)); % Normalized
        S_lookup{test}(i,1) = alpha(i);
        S_lookup{test}(i,2) = sum(H{test}(:).*mat_temp{i}(:));
        S_lookup{test}(i,3) = sum(H{test}(:).*mat_temp2{i}(:));
    end
    if test == 3
        S_hat_lookup = S_lookup;
        for j = 1:3
            nrmlz(j) = sum(H{j}(:));
            S_hat_lookup{j}(:,2) = S_lookup{j}(:,2)/sum(H{j}(:)); % Normalized
            S_hat_lookup{j}(:,3) = S_lookup{j}(:,3)/sum(H{j}(:));
            plot(S_hat_lookup{j}(:,1),S_hat_lookup{j}(:,2))
            hold on
        end
    end
end

title('Goodness Metric Plots','interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
ylabel('$S_{metric}=\sum_{i,j}\hat{hist}\left(k_i,\lambda_j,test\right)\;.*\;\hat{S}_\alpha\left(k_i,\lambda_j,\alpha\right)$',...
    'interpreter','latex')
legend('Biaxial; $U_d = 4.5 mm$','Uniaxial; $U_d = 1 mm$','Uniaxial (MR-u); $U_d = 7 mm$',...
    'interpreter','latex')
saveas(gcf,['23-0626-Biaxial_Uniaxial_Comp/goodness_' type{ii} '.png'])
saveas(gcf,['23-0626-Biaxial_Uniaxial_Comp/goodness_' type{ii} '.pdf'])
    
end 