clc;clear
addpath(genpath('Functions'))
camposit = [90 180 105];

load('Raw Data\Comp_Mag_hCI_flipped.mat')
load('Experimental\No_Holes\refpositions.mat')
load('Experimental\No_Holes\MRI-3Ddefs_SimpleShear_220210_0020noholes.mat')
new_mask = erode_mask(mask);
%% Real component

X_pos = X{2}(1:64,:,:); temp{1} = X{1}(1:64,:,:); temp{2} = X{2}(1:64,:,:) - mean([max(X_pos(:)) min(X_pos(:))]);...
    temp{3} = X{3}(1:64,:,:); temp{4} = real(hCIz(1:64,:,:)).*new_mask(1:64,:,:);
scatterColor3(temp,51,cividis,1,0.9,51,camposit,[min(real(hCIz(:))) max(real(hCIz(:)))]); ylim([-5 5])
saveas(gcf,[pwd '\5a-Complex\real_3d.png'])

figure(); slice(temp{1},temp{2},temp{3},real(hCIz(1:64,:,:)).*new_mask(1:64,:,:),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap(cividis); caxis([min(real(hCIz(:))) max(real(hCIz(:)))]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\real_slice_mask.png'])
figure(); slice(temp{1},temp{2},temp{3},real(hCIz(1:64,:,:)),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap(cividis); caxis([min(real(hCIz(:))) max(real(hCIz(:)))]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\real_slice_nomask.png'])

figure(); imagesc(real(reshape(hCIz(1:64,16,:),[64,32]))); colormap(cividis); daspect([0.1875 1 0.75]); caxis([min(real(hCIz(:))) max(real(hCIz(:)))]); axis off
saveas(gcf,[pwd '\5a-Complex\real_th_lng.png'])
figure(); imagesc(real(reshape(hCIz(1:64,:,16),[64,32]))); colormap(cividis); daspect([0.1875 0.75 1]); caxis([min(real(hCIz(:))) max(real(hCIz(:)))]); axis off
saveas(gcf,[pwd '\5a-Complex\real_th_wd.png'])
figure(); imagesc(real(reshape(hCIz(32,:,:),[32,32]))); colormap(cividis); daspect([0.75 1 0.1875]); caxis([min(real(hCIz(:))) max(real(hCIz(:)))]); axis off
saveas(gcf,[pwd '\5a-Complex\real_wd_lng.png'])

%% Imaginary component

temp{4} = imag(hCIz(1:64,:,:)).*new_mask(1:64,:,:);
scatterColor3(temp,51,cividis,1,0.9,51,camposit,[min(imag(hCIz(:))) max(imag(hCIz(:)))]); ylim([-5 5])
saveas(gcf,[pwd '\5a-Complex\imag_3d.png'])

figure(); slice(temp{1},temp{2},temp{3},imag(hCIz(1:64,:,:)).*new_mask(1:64,:,:),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap(cividis); caxis([min(imag(hCIz(:))) max(imag(hCIz(:)))]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\imag_slice_mask.png'])
figure(); slice(temp{1},temp{2},temp{3},imag(hCIz(1:64,:,:)),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap(cividis); caxis([min(imag(hCIz(:))) max(imag(hCIz(:)))]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\imag_slice_nomask.png'])

figure(); imagesc(imag(reshape(hCIz(1:64,16,:),[64,32]))); colormap(cividis); daspect([0.1875 1 0.75]); caxis([min(imag(hCIz(:))) max(imag(hCIz(:)))]); axis off
saveas(gcf,[pwd '\5a-Complex\imag_th_lng.png'])
figure(); imagesc(imag(reshape(hCIz(1:64,:,16),[64,32]))); colormap(cividis); daspect([0.1875 0.75 1]); caxis([min(imag(hCIz(:))) max(imag(hCIz(:)))]); axis off
saveas(gcf,[pwd '\5a-Complex\imag_th_wd.png'])
figure(); imagesc(imag(reshape(hCIz(32,:,:),[32,32]))); colormap(cividis); daspect([0.75 1 0.1875]); caxis([min(imag(hCIz(:))) max(imag(hCIz(:)))]); axis off
saveas(gcf,[pwd '\5a-Complex\imag_wd_lng.png'])

%% Magnitude

temp{4} = log10(hCI_mag(1:64,:,:)).*new_mask(1:64,:,:);
scatterColor3(temp,51,colormap('gray'),1,0.9,51,camposit,[2 4]); ylim([-5 5])
saveas(gcf,[pwd '\5a-Complex\mag_3d.png'])

figure(); slice(temp{1},temp{2},temp{3},log10(hCI_mag(1:64,:,:)).*new_mask(1:64,:,:),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap('gray'); caxis([2 4]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\mag_slice_mask.png'])
figure(); slice(temp{1},temp{2},temp{3},log10(hCI_mag(1:64,:,:)),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap('gray'); caxis([2 4]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\mag_slice_nomask.png'])

figure(); imagesc(log10(reshape(hCI_mag(1:64,16,:),[64,32]))); colormap('gray'); daspect([0.1875 1 0.75]); caxis([2 4]); axis off
saveas(gcf,[pwd '\5a-Complex\mag_th_lng.png'])
figure(); imagesc(log10(reshape(hCI_mag(1:64,:,16),[64,32]))); colormap('gray'); daspect([0.1875 0.75 1]); caxis([2 4]); axis off
saveas(gcf,[pwd '\5a-Complex\mag_th_wd.png'])
figure(); imagesc(log10(reshape(hCI_mag(32,:,:),[32,32]))); colormap('gray'); daspect([0.75 1 0.1875]); caxis([2 4]); axis off
saveas(gcf,[pwd '\5a-Complex\mag_wd_lng.png'])

%% Phase

temp{4} = angle(hCIz(1:64,:,:)).*new_mask(1:64,:,:);
scatterColor3(temp,51,inferno,1,0.9,51,camposit,[-pi pi]); ylim([-5 5])
saveas(gcf,[pwd '\5a-Complex\phase_3d.png'])

figure(); slice(temp{1},temp{2},temp{3},angle(hCIz(1:64,:,:)).*new_mask(1:64,:,:),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap(inferno); caxis([-pi pi]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\phase_slice_mask.png'])
figure(); slice(temp{1},temp{2},temp{3},angle(hCIz(1:64,:,:)),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap(inferno); caxis([-pi pi]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\phase_slice_nomask.png'])

figure(); imagesc(angle(reshape(hCIz(1:64,16,:),[64,32]))); colormap(inferno); daspect([0.1875 1 0.75]); caxis([-pi pi]); axis off
saveas(gcf,[pwd '\5a-Complex\phase_th_lng.png'])
figure(); imagesc(angle(reshape(hCIz(1:64,:,16),[64,32]))); colormap(inferno); daspect([0.1875 0.75 1]); caxis([-pi pi]); axis off
saveas(gcf,[pwd '\5a-Complex\phase_th_wd.png'])
figure(); imagesc(angle(reshape(hCIz(32,:,:),[32,32]))); colormap(inferno); daspect([0.75 1 0.1875]); caxis([-pi pi]); axis off
saveas(gcf,[pwd '\5a-Complex\phase_wd_lng.png'])

%% Strain
E_l = F_t_to_E_l(F_t,'3D');
for t = 1:3
    for i = 1:3
        for j = 1:3
            E_l{t}{i,j} = -E_l{t}{i,j};
        end
    end
end

temp{4} = E_l{2}{3,3}(1:64,:,:).*new_mask(1:64,:,:);
scatterColor3(temp,51,colormap(cbrewer('div','RdYlBu',51,'linear')),1,0.9,51,camposit,[-0.8 0.8]); ylim([-5 5])
saveas(gcf,[pwd '\5a-Complex\strain_3d.png'])

figure(); slice(temp{1},temp{2},temp{3},E_l{2}{3,3}(1:64,:,:).*new_mask(1:64,:,:),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap(cbrewer('div','RdYlBu',51,'linear')); caxis([-0.8 0.8]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\strain_slice_mask.png'])
figure(); slice(temp{1},temp{2},temp{3},E_l{2}{3,3}(1:64,:,:),0,-0.7,0); daspect([1 1 1]); shading interp; campos([camposit]);
colormap(cbrewer('div','RdYlBu',51,'linear')); caxis([-0.8 0.8]); axis off;
xlim([-10 10]); ylim([-5 5]); zlim([-15 15]);
saveas(gcf,[pwd '\5a-Complex\strain_slice_nomask.png'])

figure(); imagesc(reshape(E_l{2}{3,3}(1:64,16,:),[64,32])); colormap(cbrewer('div','RdYlBu',51,'linear')); daspect([0.1875 1 0.75]); caxis([-0.8 0.8]); axis off
saveas(gcf,[pwd '\5a-Complex\strain_th_lng.png'])
figure(); imagesc(reshape(E_l{2}{3,3}(1:64,:,16),[64,32])); colormap(cbrewer('div','RdYlBu',51,'linear')); daspect([0.1875 0.75 1]); caxis([-0.8 0.8]); axis off
saveas(gcf,[pwd '\5a-Complex\strain_th_wd.png'])
figure(); imagesc(reshape(E_l{2}{3,3}(32,:,:),[32,32])); colormap(cbrewer('div','RdYlBu',51,'linear')); daspect([0.75 1 0.1875]); caxis([-0.8 0.8]); axis off
saveas(gcf,[pwd '\5a-Complex\strain_wd_lng.png'])

close all
