%This script is the run function for the image processing related to the
%ligament stretching project.

clear all; %close all;
%gets rid of the view3dgui plots
%delete(findall(0));

%%%%User-Defined Inputs%%%%
%number of the reference image as a string
fnum = '1245';
%tag for output files, e.g. 'nofilt', 'prekfilt', 'postfilt'
tag = 'CD';
%camera position/angle for the 3d plots
cp = [90 180 105];
%chosen z value for plots
z = 48; %48 for ligaments, 63 for dogbone

%toggle on/off options:
%save unwrap files, mechanical fields
saving = false;
%view 3d plots
pl3d = false;
%filter before unwrapping
prefilt = true;
hfilt = [3 3 3];%[1.2*1.5 3 3];%[3 4 5];
%Mask threshold on a log scale
maskthresh = -0.6;
reps = 0;
%filter u after unwrapping
ufilt = false;
%use weighted reliability algorithm, which fills in edges last
weighting = true;
%only include the largest single connected mask object
largestConnectedObjectMasking=true;
%open view3dgui plots
plView3Ds = false;
%%%%%%%%

hstr = [num2str(hfilt(1)) num2str(hfilt(2)) num2str(hfilt(3))];

%loads the MR info
fn = glob(['DESTE*' fnum '*.mat']);
load(fn{1});
L = lambda;
try
    hsc = [abs(HIRES.axis1(2)-HIRES.axis1(1)), abs(HIRES.axis2(2)-HIRES.axis2(1)), abs(HIRES.axis3(2)-HIRES.axis3(1))]
catch
    hfn = glob(['SHIRES*' fnum '*.mat']);
    try
        load(hfn{1});
        HIRES = struct();
        HIRES.axis1 = axis1; HIRES.axis2 = axis2; HIRES.axis3 = axis3; HIRES.magnitude = magnitude;
    catch
        HIRES = struct();
        HIRES.axis1 = axis1; HIRES.axis2 = axis2; HIRES.axis3 = axis3; HIRES.magnitude = abs(data);
    end
    hsc = [abs(HIRES.axis1(2)-HIRES.axis1(1)), abs(HIRES.axis2(2)-HIRES.axis2(1)), abs(HIRES.axis3(2)-HIRES.axis3(1))]
end
%reload to overwrite any vars with the same name
load(fn{1});

%Make a mask out of the HIRES info
hires_mag = HIRES.magnitude;%flip(HIRES.magnitude,3);
[hx,hy,hz] = meshgrid(HIRES.axis2,HIRES.axis1,HIRES.axis3);
try
    CI_mag = (abs(data(:,:,:,1))+abs(data(:,:,:,2))+abs(data(:,:,:,3)))/3;
catch
    data = cat(4,data,data,data);
    CI_mag = (abs(data(:,:,:,1))+abs(data(:,:,:,2))+abs(data(:,:,:,3)))/3;
end
[cx,cy,cz] = meshgrid(axis2,axis1,axis3);

hCI_mag = interp3(cx,cy,cz,CI_mag,hx,hy,hz);
mask = nan(size(hCI_mag)); %mask(hires_mag>2E4) = 1;
%start with some baseline thresholding
mask(log10(hCI_mag/max(hCI_mag(:)))>maskthresh) = 1; %e.g. -0.6 for Dragon Skin, -1.4 for ligs
%mask(log10(hires_mag/max(hires_mag(:)))>-0.6) = 1; %-0.6 for Dragon Skin, -1.4 for ligs
if plView3Ds, view3dgui(mask); end
if largestConnectedObjectMasking
    %view3dgui(mask);
    r = regionprops(mask>0);
    bwc = bwconncomp(mask>0);
    %find the largest connected object...
    for k=1:length(r)
        area(k) = r(k).Area;
    end
    [MA, idx] = max(area);
    mask = mask - mask; %makes it zeros again;
    %then remove all the non-connected smaller objects
    mask(bwc.PixelIdxList{idx}) = 1; mask(mask==0) = nan;
    try 
    %    load('mask.mat');
    %catch
        load('mask_HIRES.mat');
        mask = round(mask_HIRES(:,1:2:127,:)/2+mask_HIRES(:,2:2:128,:)/2);
    end
end

%Look at the mask if plot view 3D is on
if plView3Ds
    view3dgui(mask);
end
%mask(hCI_mag>0) = 1;


if prefilt
    
    %hCIx = ndnanfilter(hCIx,'hamming',hfilt);
    %hCIy = ndnanfilter(hCIy,'hamming',hfilt);
    %hCIz = ndnanfilter(hCIz,'hamming',hfilt);
    
    datax = blur3d(data(:,:,:,1),'vox',hfilt);
    datay = blur3d(data(:,:,:,2),'vox',hfilt);
    dataz = blur3d(data(:,:,:,3),'vox',hfilt);
    
    if reps ~=0
        for m=1:reps
            datax = blur3d(datax,'vox',hfilt);
            datay = blur3d(datay,'vox',hfilt);
            dataz = blur3d(dataz,'vox',hfilt);
        end
    end
    
end

hCIx = interp3(cx,cy,cz,datax,hx,hy,hz);
hCIy = interp3(cx,cy,cz,datay,hx,hy,hz);
hCIz = interp3(cx,cy,cz,dataz,hx,hy,hz);

chy = 20;
figure(31);
subplot(1,4,1);
imagesc(squeeze(real(hCIx(:,chy,:)))); colormap cividis; axis image;
set(gca,'Ydir','normal');
subplot(1,4,2);
imagesc(squeeze(imag(hCIx(:,chy,:)))); colormap cividis; axis image;
set(gca,'Ydir','normal');
subplot(1,4,3);
imagesc(squeeze(abs(hCIx(:,chy,:)))); colormap cividis; axis image;
set(gca,'Ydir','normal');
subplot(1,4,4);
imagesc(squeeze(angle(hCIx(:,chy,:).*mask(:,chy,:)))); colormap cividis; axis image;
set(gca,'Ydir','normal');

shp = alphaShape(hx(mask==1),hy(mask==1),hz(mask==1))
%figure; plot(shp);

%Fix the mask to remove internal holes using the alphashape
TF = inShape(shp,hx,hy,hz); mask_ = nan(size(TF)); mask_(TF) = 1;
mask = mask_;


%% Check if there's a separate mask

%  try load('vol_mask.mat')%load([fnum '_ref_vol_mask.mat'])
%      mask = nan(size(tf));
%      mask(tf) = 1;
%      shp = alphaShape(hx(mask==1),hy(mask==1),hz(mask==1))
%  end

%Signal histogram from the MR - blue = all points, orange = those in mask
figure; histogram(log10(hCI_mag/max(hCI_mag(:))),-3:0.01:0)
hold on; histogram(log10(mask.*hCI_mag/max(hCI_mag(:))),-3:0.01:0);


%Temporarily removed this for strain shift calcs, 180810
%hCIx = hCIx.*mask; hCIy = hCIy.*mask; hCIz = hCIz.*mask;
HCI = cat(4,hCI_mag/max(hCI_mag(:)),hCI_mag/max(hCI_mag(:)),hires_mag/max(hires_mag(:)));

%Complex plane plot of the phoase info - large radius = more signal
figure(111); scatter(real(hCIx(mask==1)), imag(hCIx(mask==1)),'.'); axis image;
%subplot(1,2,2); scatter(real(hCIx(mask==1)), imag(hCIx(mask==1)),'.'); axis image;

%% Error estimation

mask_xlims = [80 110]; 
invmask = isnan(mask); invmask(1:mask_xlims(1),:,:) = false; invmask(mask_xlims(2):end,:,:) = false;


%xyz_complexError = [std(real(datax(invmask))), std(real(datay(invmask))), std(real(dataz(invmask)));...
%                    std(imag(datax(invmask))), std(imag(datay(invmask))), std(imag(dataz(invmask)))];
%err_circ = sqrt(mean(xyz_complexError(1,:))^2+mean(xyz_complexError(2,:))^2);
err_circ = 10.^nanmedian(log10(abs(hCI_mag(invmask))));

rdatax = abs(datax); rdatay = abs(datay); rdataz = abs(dataz);
err_Lx = sqrt(rdatax.^2-err_circ.^2); err_Ly = sqrt(rdatay.^2-err_circ.^2); err_Lz = sqrt(rdataz.^2-err_circ.^2);
err_Lx(imag(err_Lx)>0) = nan; err_Lx = abs(err_Lx);
err_Ly(imag(err_Ly)>0) = nan; err_Ly = abs(err_Ly);
err_Lz(imag(err_Lz)>0) = nan; err_Lz = abs(err_Lz);

datax_err = 2*atan(err_circ/err_Lx); datax_err(isnan(datax_err)) = pi;
datay_err = 2*atan(err_circ/err_Ly); datay_err(isnan(datay_err)) = pi;
dataz_err = 2*atan(err_circ/err_Lz); dataz_err(isnan(dataz_err)) = pi;

%% Perform the complex divide filtering

LResMask = interp3(hx,hy,hz,mask,cx,cy,cz);
filts{1} = 'optimal3'; %filts{2} = 'optimal3'; %filts{3} = 'optimal5';

for j = 1:length(filts)
    
    filt = filts{j};
    [Qij{1,1}(:,:,:), Qij{1,2}(:,:,:), Qij{1,3}(:,:,:)] = complexDivideN(hCIx.*mask,filt);
    [Qij{2,1}(:,:,:), Qij{2,2}(:,:,:), Qij{2,3}(:,:,:)] = complexDivideN(hCIy.*mask,filt);
    [Qij{3,1}(:,:,:), Qij{3,2}(:,:,:), Qij{3,3}(:,:,:)] = complexDivideN(hCIz.*mask,filt);
    
    % [oQij{1,1}(:,:,:), oQij{1,2}(:,:,:), oQij{1,3}(:,:,:)] = complexDivideN(data(:,:,:,1),filt);
    % [oQij{2,1}(:,:,:), oQij{2,2}(:,:,:), oQij{2,3}(:,:,:)] = complexDivideN(data(:,:,:,2),filt);
    % [oQij{3,1}(:,:,:), oQij{3,2}(:,:,:), oQij{3,3}(:,:,:)] = complexDivideN(data(:,:,:,3),filt);
    
    [oQij{1,1}(:,:,:), oQij{1,2}(:,:,:), oQij{1,3}(:,:,:)] = complexDivideN(data(:,:,:,1).*LResMask,filt);
    [oQij{2,1}(:,:,:), oQij{2,2}(:,:,:), oQij{2,3}(:,:,:)] = complexDivideN(data(:,:,:,2).*LResMask,filt);
    [oQij{3,1}(:,:,:), oQij{3,2}(:,:,:), oQij{3,3}(:,:,:)] = complexDivideN(data(:,:,:,3).*LResMask,filt);
    
    delta(1) = axis1(2)-axis1(1);
    delta(2) = axis2(2)-axis2(1);
    delta(3) = axis3(2)-axis3(1);
    
    %Create GUij == u_i,j
    GUij = cell(3,3);
    for e =1:3
        for g =1:3
            try
                GUij{e,g}(:,:,:) = -Qij{e,g}(:,:,:)/(2*pi)*lambda(e)/(2*delta(g));
            catch
                GUij{e,g}(:,:,:) = -Qij{e,g}(:,:,:)/(2*pi)*lambda/(2*delta(g));
            end
        end
    end
    
%     E11_err = -datax_err/(2*pi)*lambda(1)/(2*delta(1));
%     E12_err = -datax_err/(2*pi)*lambda(1)/(2*delta(2));
%     E22_err = -datay_err/(2*pi)*lambda(2)/(2*delta(2));
%     E23_err = -datax_err/(2*pi)*lambda(2)/(2*delta(3));
%     E13_err = -datax_err/(2*pi)*lambda(1)/(2*delta(3));
%     E33_err = -dataz_err/(2*pi)*lambda(3)/(2*delta(3));
    
%    figure (91); 
    %subplot(3,3,1);
%    histogram(E11_err(:).*mask(:)); hold on; histogram(Eij{1,1}(:)*mask(:))
%     subplot(3,3,2);
%     histogram(E12_err(:).*mask(:)); hold on; histogram(Eij{1,2}(:)*mask(:))
%     subplot(3,3,2);
%     histogram(E12_err(:).*mask(:)); hold on; histogram(Eij{1,3}(:)*mask(:))

    
for i=1:3
    for j=1:3
        Fij{i,j} = GUij{i,j}+dij(i,j);
    end
end

    %Calculate small strain eps, Lagrange-Green strain Eij, and effective strain Em
    epsij = calculateEpsij(GUij);
    Eij = calculateEij(GUij);
    Em = calculateEm(Eij);
    
    if ~prefilt
        hstr = 'off';
    end
    
    %Save into a file
    save(['Fij_' filt 'h' hstr '_' fnum '.mat'], 'epsij', 'Eij', 'Fij', 'Em','mask', 'hfilt');
end

%Strain histogram
figure(80); hold on; histogram(Eij{1,1})
histogram(Eij{2,2}); histogram(Eij{3,3});

%% Strain plots in 3D using my function scatterColor3

if ~pl3d
    %magmaVals = colormap(magma);
    rbVals = flipud(brewermap(256,'RdYlBu'));
    XYZC{1} = hx; XYZC{2} = hy; XYZC{3} = hz;
    XYZC{4} = Eij{1,1}.*mask;
    %XYZC{4} = smooth3(Eij{1,1},'gaussian',[5 3 3]).*mask;%Eij{ij(1),ij(2)}.*mask;% - (max(-ui{comp}(:).*mask(:))-min(-ui{comp}(:).*mask(:)))/2;
    %XYZC{4} = ndnanfilter(Em,'hamming',[5 3 3]).*mask;
    
    fxx = scatterColor3(XYZC, 21, rbVals, 1.0, 0.6, 20, cp, [-0.2 0.2]);
    %hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.01; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
    %hold off;
    %zlim([0 10])
    %xlim([-5 5])
    %ylim([-25 15])
    %saveas(gcf,'wiggleE11.png')
    
    XYZC{4} = Eij{2,2}.*mask;
    fyy = scatterColor3(XYZC, 21, rbVals, 1.0, 0.6, 20, cp, [-0.2 0.2]);
    %hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.01; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
    %hold off;
    %zlim([0 10])
    %xlim([-5 5])
    %ylim([-25 15])
    %saveas(gcf,'wiggleE22.png')
    
    XYZC{4} = Eij{3,3}.*mask;
    fzz = scatterColor3(XYZC, 21, rbVals, 1.0, 0.6, 20, cp, [-0.2 0.2]);
    %hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.01; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
    %hold off;
    %zlim([0 10])
    %xlim([-5 5])
    %ylim([-25 15])
    %saveas(gcf,'wiggleE33.png')
    
    rbVals = colormap(parula);
    XYZC{4} = Em.*mask;
    fm = scatterColor3(XYZC, 21, rbVals, 1.0, 0.6, 12, cp, [0 0.4]);
    %hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.01; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
    %hold off;
    %zlim([0 10])
    %xlim([-5 5])
    %ylim([-25 15])
    %saveas(gcf,'wiggleEm.png')
end

%% Temporary code for analyzing the QC stuff

% QCamp = abs(QC.data(:,:,:,1,1))./abs(QC.data(:,:,:,2,1));
% QCin = QCamp(:,1:2:63,1:2:63,1,1).*mask(:,1:2:63,1:2:63);
% view3dgui(QCin)

%% Calculating displacements using unwrapping algorithm

if weighting
    w = ones(size(hCIx));
    w = w.*(hCI_mag/max(hCI_mag(:))).^2;
else
    w = ones(size(hCIx));
end
%try weighting by the phase information alone first

tic
ui{1} = WeightedUnwrap3(angle(hCIx).*mask,w,mask)*lambda(1)/(2*pi);
toc
tic
ui{2} = WeightedUnwrap3(angle(hCIy).*mask,w,mask)*lambda(2)/(2*pi);
toc
tic
ui{3} = WeightedUnwrap3(angle(hCIz).*mask,w,mask)*lambda(3)/(2*pi);
toc

if ~pl3d
    view3dgui(ui{1})
    view3dgui(ui{2})
    view3dgui(ui{3})
end

%% Calculate the invariants to plot on I1-I2 plot



for i=1:3
    for j=1:3
        Bij{i,j} = Fij{i,1}.*Fij{j,1} + Fij{i,2}.*Fij{j,2} + Fij{i,3}.*Fij{j,3};
    end
end

Bii = Bij{1,1} + Bij{2,2} + Bij{3,3};
BijBji = Bij{1,1}.*Bij{1,1} + Bij{1,2}.*Bij{2,1} + Bij{1,3}.*Bij{3,1}...
    + Bij{2,1}.*Bij{1,2} + Bij{2,2}.*Bij{2,2} + Bij{2,3}.*Bij{3,2}...
    + Bij{3,1}.*Bij{1,3} + Bij{3,2}.*Bij{2,3} + Bij{3,3}.*Bij{3,3};

I1 = Bii;
I2 = 0.5*(Bii.^2 - BijBji);
J = calculateJ(Fij);

I1_ = I1.*J.^(-2/3);
I2_ = I2.*J.^(-4/3);
%I2__ = Bij{1,1}.*Bij{2,2}+Bij{2,2}.*Bij{3,3}+Bij{1,1}.*Bij{3,3}...
%    -(Bij{1,2}.*Bij{2,1}+Bij{1,3}.*Bij{3,1}+Bij{2,3}.*Bij{3,2});

figure(30);
plot(I1_(:),I2_(:),'.');
xlim([3 3.1]);
ylim([3 3.1]);
hold on;

lam = 0.5:0.01:3;
I1lam = lam.^2+2./lam;
I1uniax_lam = lam.^2+2./lam;
I2uniax_lam = 2*lam+(lam.^-2);
I1ps_lam = 1+lam.^2+(lam.^-2);
I2ps_lam = I1ps_lam;
I1biax_lam = 2*lam.^2+lam.^-4;
I2biax_lam = lam.^4+2./lam.^2;
plot(I1uniax_lam,I2uniax_lam,'-','color','r');
plot(I1ps_lam,I2ps_lam,'-','color','r');
plot(I1biax_lam,I2biax_lam,'-','color','r');




%% Old: Plot comparison distributions for different lambdas
% str = '1248';
% magmavals = colormap(magma);
%
% figure;
% for p=1:4
%     files = glob(['G:\My Drive\RESEARCH\Ligament Stretcher\Data\dragonskin_nicked_082418\dragonskin_nicked01\nicked_lambda' str(p) '_stretch1\*.mat']);
%     load(files{1}, 'strains');
%     u11 = strains(:,:,:,1,1).*mask(:,:,:);
%     hold on; [k, xi] = ksdensity(u11(:));
%     plot(xi,k,'LineWidth',2,'Color',magmavals(1+floor((4-p)/4*64),:));
%
% end
% title('Direct Complex Division');
