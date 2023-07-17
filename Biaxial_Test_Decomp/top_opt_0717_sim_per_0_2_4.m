function output = top_opt_0717_sim_per_0_2_4()
% Decomposes an input file for biaxial test into K2/K3 data.
% Note that all length parameters are in mm.
warning off
addpath(genpath('Data'))
addpath(genpath('Functions'))

parallel = false; % parfor loop for the decomposition (timing is the same)
bin_res = 0.01; % Bin resolution for histograms
maxLam = 2; % Maximum limit for lambda in plots
clim = [-5 -1.5]; % Colorbar limits for galaxy plots

% Def grad storage
for i = 1:3
    if i == 1
        load('Data/23-0717-sim-Sinwaves/Eco_sq-8mm_sin-per-0_sin-amp-0mm_tet.mat');
        [k3_temp,k2_temp,~] = param_decoup_main(output.F_t,parallel,'ftok2andk3');
    elseif i == 2
        load('Data/23-0717-sim-Sinwaves/Eco_sq-8mm_sin-per-2_sin-amp-1mm_tet.mat');
        [k3_temp,k2_temp,~] = param_decoup_main(output.F_t,parallel,'ftok2andk3');
    elseif i == 3
        load('Data/23-0717-sim-Sinwaves/Eco_sq-8mm_sin-per-4_sin-amp-1mm_tet.mat');
        [k3_temp,k2_temp,~] = param_decoup_main(output.F_t,parallel,'ftok2andk3');
    end
    k3{i} = k3_temp{1};
    k2{i} = k2_temp{1};
end

output.k3 = k3; output.k2 = k2;

output.lbl{1} = 'd=8mm;N=0;A=0mm';
output.lbl{2} = 'd=8mm;N=2;A=1mm';
output.lbl{3} = 'd=8mm;N=4;A=1mm';

mu = 1; % Normalizing mu
alpha = linspace(-5,5,1001);

k_2 = 0:bin_res:floor(log(maxLam)/sqrt(2/3)/bin_res)*bin_res;
k_3=-1:2*bin_res:1;

g1 = @(k_3) sqrt(2/3)*sin(-asin(k_3)/3+(2*pi/3));
g2 = @(k_3) sqrt(2/3)*sin(-asin(k_3)/3);
g3 = @(k_3) sqrt(2/3)*sin(-asin(k_3)/3-(2*pi/3));

Smu_gen = @(k_3) g1(k_3)*exp(g1(k_3)*bsxfun(@times,alpha,k_2')) + ...
    g2(k_3)*exp(g2(k_3)*bsxfun(@times,alpha,k_2')) + ...
    g3(k_3)*exp(g3(k_3)*bsxfun(@times,alpha,k_2'));
Salpha_gen = @(k_3) mu*((g1(k_3)^2)*bsxfun(@times,k_2',exp(g1(k_3)*...
    bsxfun(@times,alpha,k_2'))) + (g2(k_3)^2)*bsxfun(@times,k_2',...
    exp(g2(k_3)*bsxfun(@times,alpha,k_2'))) + (g3(k_3)^2)*bsxfun(@times,...
    k_2',exp(g3(k_3)*bsxfun(@times,alpha,k_2'))));

for i = 1:length(k_3)
    Salpha_all(:,:,i) = Salpha_gen(k_3(i));
    Smu_all(:,:,i) = Smu_gen(k_3(i));
end
output.Salpha_all = Salpha_all;
output.Smu_all = Smu_all;
clear H

for i = 1:length(output.k3)
    figure()
    [~,~,h] = ndhist(output.k2{i}(:),output.k3{i}(:),'axis',[min(k_2) max(k_2) min(k_3) max(k_3)],'fixed_bin_res',bin_res,'scale_y',2);
    temp_H{i} = h';
    imagesc([min(k_2) max(k_2)],[min(k_3) max(k_3)],log10(h/sum(h(:))));
    axis square
    colormap turbo;
    set(gca,'YDir','normal'); 
    ylabel('$k_3$, Strain state','interpreter','latex')
    xlabel('$k_2$, Log strain amplitude','interpreter','latex')
    title(['k2 vs k3 Decomposition $' output.lbl{i} '$'],'interpreter','latex')
    c = colorbar; caxis(clim);
    c.Label.String = 'log(h/sum(h(:)))';
    if i == 1
        saveas(gcf,'Data/23-0717-sim-Sinwaves/2Dhist_sq8_per0_amp0_k2_k3.png')
        saveas(gcf,'Data/23-0717-sim-Sinwaves/2Dhist_sq8_per0_amp0_k2_k3.pdf')
    elseif i == 2
        saveas(gcf,'Data/23-0717-sim-Sinwaves/2Dhist_sq8_per2_amp1_k2_k3.png')
        saveas(gcf,'Data/23-0717-sim-Sinwaves/2Dhist_sq8_per2_amp1_k2_k3.pdf')
    elseif i == 3
        saveas(gcf,'Data/23-0717-sim-Sinwaves/2Dhist_sq8_per4_amp1_k2_k3.png')
        saveas(gcf,'Data/23-0717-sim-Sinwaves/2Dhist_sq8_per4_amp1_k2_k3.pdf')
    end
end

output.H = temp_H;
H = output.H; k2 = output.k2; k3 = output.k3;

save('Data/23-0717-sim-Sinwaves/sensitivity.mat','H','k2','k3')

figure();
for test = 1:length(output.k3)
    for i = 1:length(alpha)
        mat_temp{i} = squeeze(Salpha_all(:,i,:));
        mat_temp2{i} = squeeze(Smu_all(:,i,:));
        % mat_temp{i} = mat_temp{i}./sum(mat_temp{i}(:)); % Normalized
        S_lookup{test}(i,1) = alpha(i);
        S_lookup{test}(i,2) = sum(output.H{test}(:).*mat_temp{i}(:));
        S_lookup{test}(i,3) = sum(output.H{test}(:).*mat_temp2{i}(:));
    end
    if test == length(output.k3)
        S_hat_lookup = S_lookup;
        for j = 1:length(output.k3)
            nrmlz(j) = sum(output.H{j}(:));
            S_hat_lookup{j}(:,2) = S_lookup{j}(:,2)/sum(output.H{j}(:)); % Normalized
            S_hat_lookup{j}(:,3) = S_lookup{j}(:,3)/sum(output.H{j}(:)); % Normalized
            plot(S_hat_lookup{j}(:,1),S_hat_lookup{j}(:,2),'Color',[j/length(output.k3) 0 0],'LineWidth',2)
            hold on
        end
    end
end
title('$S_{metric}$(test,$\alpha$)','interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
ylabel('$S_{metric}=\sum_{i,j}\hat{hist}\left(k_i,\lambda_j,test\right)\;.*\;\hat{S}_\alpha\left(k_i,\lambda_j,\alpha\right)$',...
    'interpreter','latex')
legend(output.lbl{1},output.lbl{2},output.lbl{3},'interpreter','latex')
saveas(gcf,'Data/23-0717-sim-Sinwaves/goodness_curve.png')
saveas(gcf,'Data/23-0717-sim-Sinwaves/goodness_curve.pdf')