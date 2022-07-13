clc;clear
addpath(genpath('Functions'))
camposit = [90 180 105];
timstp = 2;

load('Experimental\No_Holes\MRI-3Ddefs_SimpleShear_220210_0020noholes.mat')
load('Experimental\No_Holes\refpositions.mat')
new_mask = erode_mask(mask);

%% Slice of U_3

for t = 1:3
    for i = 1:3
        U_t{t}{i} = -U_t{t}{i};
    end
end
X_pos_1 = X{1}(1:64,:,:).*mask(1:64,:,:); 
X_pos_2 = X{2}(1:64,:,:).*mask(1:64,:,:); 
X_pos_3 = X{3}(1:64,:,:).*mask(1:64,:,:);
temp{1} = X{1}(1:64,:,:) - mean([max(X_pos_1(:)) min(X_pos_1(:))]);
temp{2} = X{2}(1:64,:,:) - mean([max(X_pos_2(:)) min(X_pos_2(:))]);
temp{3} = X{3}(1:64,:,:) - mean([max(X_pos_3(:)) min(X_pos_3(:))]);
temp{4} = U_t{timstp}{3}(1:64,:,:).*new_mask(1:64,:,:);

figure(); slice(temp{1},temp{2},temp{3},temp{4},0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap(cbrewer('div','RdYlBu',51,'linear')); caxis([-prescribedU(timstp) prescribedU(timstp)]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5bcd_V2-Slice_U\slice.png'])

%% Experimental w/ No Holes

% U1
temp{4} = U_t{timstp}{1}(1:64,:,:).*mask(1:64,:,:);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Exp_NoHoles_U1.png']);
% U2
temp{4} = U_t{timstp}{2}(1:64,:,:).*mask(1:64,:,:);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Exp_NoHoles_U2.png']);
% U3
temp{4} = U_t{timstp}{3}(1:64,:,:).*mask(1:64,:,:);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Exp_NoHoles_U3.png']);

%% Experimental w/ Holes

load('Experimental\Holes\MRI-3Ddefs_SimpleShear_220307_0020_5MMholes_1300Tol.mat')
load('Experimental\Holes\refpositions.mat')
new_mask = erode_mask(mask);

for t = 1:2
    for i = 1:3
        U_t{t}{i} = -U_t{t}{i}; % Changes RdYlBu colors
    end
end
X_pos_1 = X{1}(1:64,:,:).*mask(1:64,:,:); 
X_pos_2 = X{2}(1:64,:,:).*mask(1:64,:,:); 
X_pos_3 = -X{3}(1:64,:,:).*mask(1:64,:,:); % Reorients sample
temp{1} = X{1}(1:64,:,:) - mean([max(X_pos_1(:)) min(X_pos_1(:))]);
temp{2} = X{2}(1:64,:,:) - mean([max(X_pos_2(:)) min(X_pos_2(:))]);
temp{3} = -X{3}(1:64,:,:) - mean([max(X_pos_3(:)) min(X_pos_3(:))]);
temp{4} = U_t{timstp}{3}(1:64,:,:).*mask(1:64,:,:);

% U1
temp{4} = U_t{timstp}{1}(1:64,:,:).*mask(1:64,:,:);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Exp_Holes_U1.png']);
% U2
temp{4} = U_t{timstp}{2}(1:64,:,:).*mask(1:64,:,:);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Exp_Holes_U2.png']);
% U3
temp{4} = U_t{timstp}{3}(1:64,:,:).*mask(1:64,:,:);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Exp_Holes_U3.png']);

%% Simulation w/ NoHoles
load('Simulations\No_Holes\MRI-3Ddefs_SimpleShear_220303_SimulationNoHoles.mat')
load('Simulations\No_Holes\refpositions.mat')

for t = 1:3
    for i = 1:3
        U_t{t}{i} = -U_t{t}{i}; % Changes RdYlBu colors
    end
end

temp{1} = reshape(X_ref{2}-0.5*max(X_ref{2}),[17,25,50]); temp{2} = -reshape(X_ref{1}-0.5*max(X_ref{1}),[17,25,50]);...
    temp{3} = reshape(X_ref{3}-0.5*max(X_ref{3}),[17,25,50]);
% U1
temp{4} = reshape(U_t{timstp}{1},[17,25,50]);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Sim_NoHoles_U1.png']);
% U2
temp{4} = reshape(U_t{timstp}{2},[17,25,50]);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Sim_NoHoles_U2.png']);
% U3
temp{4} = reshape(U_t{timstp}{3},[17,25,50]);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Sim_NoHoles_U3.png']);

%% Simulation w/ Holes
load('Simulations\Holes\MRI-3Ddefs_SimpleShear_220305_SimulationHoles.mat')
load('Simulations\Holes\refpositions.mat')

for t = 1:2
    for i = 1:3
        U_t{t}{i} = -U_t{t}{i};
    end
end

temp{1} = reshape(X_ref{2}-0.5*max(X_ref{2}),[2,17,1319]); temp{2} = -reshape(X_ref{1}-0.5*max(X_ref{1}),[2,17,1319]);...
    temp{3} = reshape(X_ref{3}-0.5*max(X_ref{3}),[2,17,1319]);
% U1
temp{4} = reshape(U_t{timstp}{1},[2,17,1319]);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Sim_Holes_U1.png']);
% U2
temp{4} = reshape(U_t{timstp}{2},[2,17,1319]);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Sim_Holes_U2.png']);
% U3
temp{4} = reshape(U_t{timstp}{3},[2,17,1319]);...
    scatterColor3(temp,51,cbrewer('div','RdYlBu',51,'linear'),1,0.8,50,camposit,[-prescribedU(timstp) prescribedU(timstp)]); ylim([-5 5])
saveas(gcf,[pwd '\5bcd_V2-Slice_U\Sim_Holes_U3.png']);

close all