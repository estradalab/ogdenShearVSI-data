%Image processing methods for ACL

clear all; close all;
cd('G:\My Drive\RESEARCH\Ligament Stretcher\Data\20180615\MAT');

f = glob('*1344*.mat');
load(f{1});

Im = magnitude;%abs(t2star);
%view3dgui(log(magnitude))
sizeIm = size(Im);
selx = 1:sizeIm(1); %[80:115 125:160];
sely = 1:sizeIm(2); %32:60;
selz = 1:sizeIm(3); %52:96;

x = 128;
y = 42;
z = 72;

pick = 'z';

Ix = squeeze(Im(x,sely,selz));
Iy = squeeze(Im(selx,y,selz));
Iz = squeeze(Im(selx,sely,z));

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
% figure;
imagesc(I); colormap magma; axis image; %caxis([0 0.03]);
set(gca,'YDir','normal');
BW = edge(I,'Canny',0.05);

% figure(8); hold on; histogram(t2star(selx,sely,selz),0:0.001:0.1)
% selt2star = t2star(selx,sely,selz);
% median(selt2star(:))
%%

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

fnRO = glob('*1120*');
load(fnRO{1});
RO = compleximage;
L(1) = lambda*10^-3;

fnPE = glob('*1132*');
load(fnPE{1});
PE = compleximage;
L(2) = lambda*10^-3;

fnSL = glob('*1143*');
load(fnSL{1});
SL = compleximage;
L(3) = lambda*10^-3;

% switch fn
%     case 'displacement1mm'
%         data = displacement1mm;
%     case 'displacement2mm'
%         data = displacement2mm;
%     case 'displacement3mm'
%         data = displacement3mm;
% end

sc = [axis1(2)-axis1(1), axis2(2)-axis2(1), axis3(2)-axis3(1)];

lambda = L; %microns per 2pi in [x y z]
ph_xA = RO;
ph_yA = PE;
ph_zA = SL;
ph_allmag = 1/3*(abs(RO)+abs(PE)+abs(SL));
%view3dgui(ph_allmag)


switch pick
    case 'x'
        ph_x = squeeze(RO(x,:,:));
        ph_y = squeeze(PE(x,:,:));
        ph_z = squeeze(SL(x,:,:));
        ph_mag = squeeze(magnitude(x,:,:));
        ph_allmagS = squeeze(ph_allmag(x,:,:));
        plsc = [sc(2) sc(3)];
    case 'y'
        ph_x = squeeze(RO(:,y,:));
        ph_y = squeeze(PE(:,y,:));
        ph_z = squeeze(SL(:,y,:));
        ph_mag = squeeze(magnitude(:,y,:));
        ph_allmagS = squeeze(ph_allmag(:,y,:));
        plsc = [sc(1) sc(3)];
    case 'z'
        ph_x = squeeze(RO(:,:,z));
        ph_y = squeeze(PE(:,:,z));
        ph_z = squeeze(SL(:,:,z));
        ph_mag = squeeze(magnitude(:,:,z));
        ph_allmagS = squeeze(ph_allmag(:,:,z));
        plsc = [sc(1) sc(2)];
end

newBW = true(size(dI));
%newBW = imbinarize(dI,mode(Im(:)/max(Im(:)))+std(Im(:)/max(Im(:))));
figure(3);
histogram(dI(:),100);

figure(2);
subplot(1,5,1);
imagesc(dI); colormap magma; axis image;
set(gca,'YDir','normal');
subplot(1,5,2);
imagesc(log(ph_allmagS)-log(max(ph_allmagS(:)))); colormap magma; axis image;
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

%figure(40)

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

%figure(30); imshow3D(ux);
%figure(31); imshow3D(permute(ux,[1 3 2]));

Iall = ph_yA;
GoldsteinUnwrapQuasi3D_pin;
uy = u_px*lambda(2)/(2*pi);

%figure(30); imshow3D(uy,[-0.5 3.5]);
%figure(31); imshow3D(permute(uy,[1 3 2]),[-0.5 3.5]);

Iall = ph_zA;
GoldsteinUnwrapQuasi3D_pin;
uz = u_px*lambda(3)/(2*pi);

% figure(30); imshow3D(uz,[-0.5 3.5]);
% figure(31); imshow3D(permute(uz,[1 3 2]),[-0.5 3.5]);


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

anag_IandPh = cat(3,dI/max(dI(:)),ph_allmagS/max(ph_allmagS(:)),ph_allmagS/max(ph_allmagS(:)));
figure(130); image(anag_IandPh)

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

