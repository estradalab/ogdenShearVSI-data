clc;clear;

%% Circular holes
% load('22-0526-Ogden_Methodical_Size\sensitivity_plts.mat')
load('Temp_Data\S_hat_lookup_VolCorr_Size_perElement.mat')
for i = 1:11
    S_Hat_alpha_circle{i} = [];
end
for i = 1:4
    for j = 1:11
        S_Hat_alpha_circle{j} = [S_Hat_alpha_circle{j} S_hat_lookup{i}((j-1)*100+1,2)];
    end
end
% load('22-0516-Ogden_Methodical\sensitivity_plts.mat')
load('Temp_Data\S_hat_lookup_VolCorr_perElement.mat')
for i = 1:11
    S_Hat_alpha_circle{i} = [S_Hat_alpha_circle{i} S_hat_lookup{2}((i-1)*100+1,2)];
end
load('Temp_Data\S_hat_lookup_VolCorr_Size_perElement.mat')
for i = 5:9
    for j = 1:11
        S_Hat_alpha_circle{j} = [S_Hat_alpha_circle{j} S_hat_lookup{i}((j-1)*100+1,2)];
    end
end

%% Square holes
load('Temp_Data\S_hat_lookup_VolCorr_Size_perElement.mat')
for i = 1:11
    S_Hat_alpha_sq{i} = [];
end
for i = [15 10 16 11 17 12 18 13 19 14]
    for j = 1:11
        S_Hat_alpha_sq{j} = [S_Hat_alpha_sq{j} S_hat_lookup{i}((j-1)*100+1,2)];
    end
end

%% Plots
indx = [1;2;3;4;5;6;7;8;9;10];

figure('Name','Sensitivity vs Size','NumberTitle','off',...
    'units','normalized','outerposition',[0.2 0.2 0.8 0.8])
for i = [1:5 7:length(S_Hat_alpha_circle)]
    if i <= 5
        subplot(2,5,i)
    else
        subplot(2,5,i-1)
    end
    plot(indx,S_Hat_alpha_circle{i}','LineWidth',2,'Color',[0.4940 0.1840 0.5560])
    xlabel('Diameter/area equivalent [mm]','interpreter','latex')
    ylabel('$$\hat{S}_\alpha$$','interpreter','latex')
    title(['$$\alpha=$$' num2str(S_hat_lookup{1}((i-1)*100+1,1))],'interpreter','latex')
    hold on
end

for i = [1:5 7:length(S_Hat_alpha_sq)]
    if i <= 5
        subplot(2,5,i)
    else
        subplot(2,5,i-1)
    end
    plot(indx,S_Hat_alpha_sq{i}','LineWidth',2,'Color',[0.3010 0.7450 0.9330])
end

sgtitle('Sensitivity vs Size of {\color[rgb]{0.4940 0.1840 0.5560}circular} or {\color[rgb]{0.3010 0.7450 0.9330}square} hole')
saveas(gcf,'22-0526-Ogden_Methodical_Size\alpha_v_size.png')