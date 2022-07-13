%Image processing methods for ACL

clear all; close all;
load('outT2star.mat');

x = 128;
y = 55;
z = 48;

pick = 'y';

Im = outT2star.amplitude;
Ix = squeeze(Im(x,:,:));
Iy = squeeze(Im(:,y,:));
Iz = squeeze(Im(:,:,z));

switch pick
    case 'x'
        I = Ix;
        plane = 'yz';
    case 'y'
        I = Iy;
        plane = 'xz';
    case 'z'
        I = Iz;
        plane = 'xy';
end


subplot(1,4,1);
imagesc(I); colormap magma; axis image;
set(gca,'YDir','normal');
BW = edge(I,'Canny',0.05);



% subplot(1,3,2);
% I = histeq(outT2star.amplitude/max(outT2star.amplitude(:)));
% imagesc(I(:,y,:)); colormap magma; axis image;
% set(gca,'YDir','normal');


subplot(1,4,2);
imagesc(BW);  colormap magma; axis image;
set(gca,'YDir','normal');

subplot(1,4,3);
net = denoisingNetwork('DnCNN');
dI = denoiseImage(I/max(I(:)), net);
imagesc(dI);  colormap magma; axis image;
set(gca,'YDir','normal');

subplot(1,4,4);
dBW = edge(dI,'Canny',0.05);
%P = poly2maskfinddBW(dBW));
RP = regionprops(dBW,'ConvexImage');
imagesc(dBW);  colormap magma; axis image;
set(gca,'YDir','normal');


%%

fn = 'displacement2mm';
load(fn);
switch fn
    case 'displacement1mm'
        data = displacement1mm;
    case 'displacement2mm'
        data = displacement2mm;
    case 'displacement3mm'
        data = displacement3mm;
end

sc = [data.axis1(2)-data.axis1(1), data.axis2(2)-data.axis2(1), data.axis3(2)-data.axis3(1)];

lambda = data.lambda; %microns per 2pi in [x y z]
ph_xA = data.phase_long_direction;
ph_yA = data.phase_across_direction;
ph_zA = data.phase_through_slice_direction;
switch pick
    case 'x'
        ph_x = squeeze(data.phase_long_direction(x,:,:));
        ph_y = squeeze(data.phase_across_direction(x,:,:));
        ph_z = squeeze(data.phase_through_slice_direction(x,:,:));
        ph_mag = squeeze(data.magnitude(x,:,:));
        plsc = [sc(2) sc(3)];
    case 'y'
        ph_x = squeeze(data.phase_long_direction(:,y,:));
        ph_y = squeeze(data.phase_across_direction(:,y,:));
        ph_z = squeeze(data.phase_through_slice_direction(:,y,:));
        ph_mag = squeeze(data.magnitude(:,y,:));
        plsc = [sc(1) sc(3)];
    case 'z'
        ph_x = squeeze(data.phase_long_direction(:,:,z));
        ph_y = squeeze(data.phase_across_direction(:,:,z));
        ph_z = squeeze(data.phase_through_slice_direction(:,:,z));
        ph_mag = squeeze(data.magnitude(:,:,z));
        plsc = [sc(1) sc(2)];
end

newBW = imbinarize(dI,mode(Im(:)/max(Im(:)))+std(Im(:)/max(Im(:))));
figure(3);
histogram(dI(:),100);
%can also use the function "angle" ffs jon why u gotta be like this

figure(2);
subplot(1,5,1);
imagesc(dI); colormap magma; axis image;
set(gca,'YDir','normal');
subplot(1,5,2);
imagesc(ph_mag); colormap magma; axis image;
set(gca,'YDir','normal');
subplot(1,5,3);
imagesc(newBW.*atan2(imag(ph_x),real(ph_x))); colormap magma; axis image;
set(gca,'YDir','normal');
subplot(1,5,4);
imagesc(newBW.*atan2(imag(ph_y),real(ph_y))); colormap magma; axis image;
set(gca,'YDir','normal');
subplot(1,5,5);
imagesc(newBW.*atan2(imag(ph_z),real(ph_z))); colormap magma; axis image;
set(gca,'YDir','normal');

%% Run the Goldstein Unwrap section
%Uncomment for 2D
% IM = ph_x;
% GoldsteinUnwrap2D_r1;
% ux = im_unwrapped*lambda(1)/(2*pi);
% 
% IM = ph_y;
% GoldsteinUnwrap2D_r1;
% uy = im_unwrapped*lambda(2)/(2*pi);
% 
% IM = ph_z;
% GoldsteinUnwrap2D_r1;
% uz = im_unwrapped*lambda(3)/(2*pi);

%Uncomment for 3D
Iall = ph_xA;
tic
GoldsteinUnwrapQuasi3D_pin;
toc
ux = u_px*lambda(1)/(2*pi);

figure(30); imshow3D(ux);
figure(31); imshow3D(permute(ux,[1 3 2]));

Iall = ph_yA;
GoldsteinUnwrapQuasi3D_pin;
uy = u_px*lambda(2)/(2*pi);

figure(30); imshow3D(uy,[-0.5 3.5]);
figure(31); imshow3D(permute(uy,[1 3 2]),[-0.5 3.5]);

Iall = ph_zA;
GoldsteinUnwrapQuasi3D_pin;
uz = u_px*lambda(3)/(2*pi);

figure(30); imshow3D(uz,[-0.5 3.5]);
figure(31); imshow3D(permute(uz,[1 3 2]),[-0.5 3.5]);


nanBW = ones(size(newBW)); nanBW(~newBW) = nan;

uL = [-0.5 1.5];

figure(9999);
subplot(1,4,1);
imagesc(dI);colormap magma; axis image;
set(gca,'YDir','normal');
title('Anatomical image, cleaned with NN');
subplot(1,4,2);
imagesc(nanBW.*ux); colormap magma; axis image; caxis(uL);
set(gca,'YDir','normal');
cc = caxis;
title('u_x, unwrapped');
subplot(1,4,3);
imagesc(nanBW.*uy); colormap magma; axis image; caxis(uL);
set(gca,'YDir','normal');
%caxis(cc);
title('u_y, unwrapped');
subplot(1,4,4);
imagesc(nanBW.*uz); colormap magma; axis image; caxis(uL);
set(gca,'YDir','normal');
%caxis(cc);
title('u_z, unwrapped');


[X1, X2] = gradientN(ux,'optimal9');
[Y1, Y2] = gradientN(uy,'optimal9');
[Z1, Z2] = gradientN(uz,'optimal9');
sX1 = X1/plsc(2); sX2 = X2/plsc(1);
sY1 = Y1/plsc(2); sY2 = Y2/plsc(1);
sZ1 = Z1/plsc(2); sZ2 = Z2/plsc(1);

%%
cL = [-0.1 0.1];

figure(10000);
subplot(1,6,1);
imagesc(nanBW.*sX1); colormap magma; axis image; caxis(cL);
set(gca,'YDir','normal');
title(['Gradient, u_x_,_' plane(2)]);
subplot(1,6,2);
imagesc(nanBW.*sX2); colormap magma; axis image; caxis(cL);
set(gca,'YDir','normal');
title(['Gradient, u_x_,_' plane(1)]);
subplot(1,6,3);
imagesc(nanBW.*sY1); colormap magma; axis image; caxis(cL);
set(gca,'YDir','normal');
title(['Gradient, u_y_,_' plane(2)]);
subplot(1,6,4);
imagesc(nanBW.*sY2); colormap magma; axis image; caxis(cL);
set(gca,'YDir','normal');
title(['Gradient, u_y_,_' plane(1)]);
subplot(1,6,5);
imagesc(nanBW.*sZ1); colormap magma; axis image; caxis(cL);
set(gca,'YDir','normal');
title(['Gradient, u_z_,_' plane(2)]);
subplot(1,6,6);
imagesc(nanBW.*sZ2); colormap magma; axis image; caxis(cL);
set(gca,'YDir','normal');
title(['Gradient, u_z_,_' plane(1)]);


switch pick
    case 'x'
 
        nY2 = nanBW.*sY2;
        disp('Small strain, eps_yy');
        nanmean(reshape(nY2,[1 numel(nY2)]))
        
        nZ1 = nanBW.*sZ1;
        disp('Small strain, eps_zz');
        nanmean(reshape(nZ1,[1 numel(nZ1)]))
                
        disp('Lagrangian strain, E_yy');
        Eyy = sY2 + 0.5*(sX2.^2+sY2.^2+sZ2.^2);
        nEyy = nanBW.*Eyy;
        nanmean(reshape(nEyy,[1 numel(nEyy)]))
        
        disp('Lagrangian strain, E_zz');
        Ezz = sZ1 + 0.5*(sX1.^2+sY1.^2+sZ1.^2);
        nEzz = nanBW.*Ezz;
        nanmean(reshape(nEzz,[1 numel(nEzz)]))
        
    case 'y'
        
        nX2 = nanBW.*sX2;
        disp('Small strain, eps_xx');
        nanmean(reshape(nX2(100:150,:),[1 numel(nX2(100:150,:))]))
        
        nZ1 = nanBW.*sZ1;
        disp('Small strain, eps_zz');
        nanmean(reshape(nZ1(100:150,:),[1 numel(nZ1(100:150,:))]))
                
        disp('Lagrangian strain, E_xx');
        Exx = sX2 + 0.5*(sX2.^2+sY2.^2+sZ2.^2);
        nExx = nanBW.*Exx;
        nanmean(reshape(nExx(100:150,:),[1 numel(nExx(100:150,:))]))
        
        disp('Lagrangian strain, E_zz');
        Ezz = sZ1 + 0.5*(sX1.^2+sY1.^2+sZ1.^2);
        nEzz = nanBW.*Ezz;
        nanmean(reshape(nEzz(100:150,:),[1 numel(nEzz(100:150,:))]))
        
    case 'z'
        
        nX2 = nanBW.*sX2;
        disp('Small strain, eps_xx');
        nanmean(reshape(nX2(100:150,:),[1 numel(nX2(100:150,:))]))
        
        nY1 = nanBW.*sY1;
        disp('Small strain, eps_yy');
        nanmean(reshape(nY1(100:150,:),[1 numel(nY1(100:150,:))]))
                
        disp('Lagrangian strain, E_xx');
        Exx = sX2 + 0.5*(sX2.^2+sY2.^2+sZ2.^2);
        nExx = nanBW.*Exx;
        nanmean(reshape(nExx(100:150,:),[1 numel(nExx(100:150,:))]))
        
        disp('Lagrangian strain, E_yy');
        Eyy = sY1 + 0.5*(sX1.^2+sY1.^2+sZ1.^2);
        nEyy = nanBW.*Eyy;
        nanmean(reshape(nEyy(100:150,:),[1 numel(nEyy(100:150,:))]))
                
end

