%This script is the run function for the image processing related to the
%ligament stretching project.

clear all; close all;
%gets rid of the view3dgui plots
%delete(findall(0));

%%%%User-Defined Inputs%%%%
%number of the reference image as a string
fnum = '1123';
%tag for output files, e.g. 'nofilt', 'prekfilt', 'postfilt'
tag = 'CD';
%camera position/angle for the 3d plots
cp = [90 180 105];
%chosen z value for plots
z = 69; %48 for ligaments, 63 for dogbone

%toggle on/off options:
%save unwrap files, mechanical fields
saving = false;
%view 3d plots
pl3d = false;
%filter before unwrapping
prefilt = false;
%filter u after unwrapping
ufilt = false;
%use weighted reliability algorithm
weighting = false;
%only include the largest single connected mask object
largestConnectedObjectMasking=true;
%open view3dgui plots
plView3Ds = false;
%flag for if Callan made a custom mask
callanmask = false;
%%%%%%%%

%loads the MR info
fn = glob(['DESTE_st*' fnum '*.mat']);
load(fn{1});

L = lambda;
try
    hsc = [abs(HIRES.axis1(2)-HIRES.axis1(1)), abs(HIRES.axis2(2)-HIRES.axis2(1)), abs(HIRES.axis3(2)-HIRES.axis3(1))]
catch
    hfn = glob(['HIRES*' fnum '*.mat']);
    load(hfn{1});
    HIRES = struct();
    HIRES.axis1 = axis1; HIRES.axis2 = axis2; HIRES.axis3 = axis3; HIRES.magnitude = magnitude;
    hsc = [abs(HIRES.axis1(2)-HIRES.axis1(1)), abs(HIRES.axis2(2)-HIRES.axis2(1)), abs(HIRES.axis3(2)-HIRES.axis3(1))]
end
%reload to overwrite any vars with the same name
load(fn{1});

%make a mask out of the HIRES info
hires_mag = HIRES.magnitude;%flip(HIRES.magnitude,3);
[hx,hy,hz] = meshgrid(HIRES.axis2,HIRES.axis1,HIRES.axis3);
CI_mag = (abs(data(:,:,:,1))+abs(data(:,:,:,2))+abs(data(:,:,:,3)))/3;
[cx,cy,cz] = meshgrid(axis2,axis1,axis3);
hCI_mag = interp3(cx,cy,cz,CI_mag,hx,hy,hz);
mask = nan(size(hires_mag)); %mask(hires_mag>2E4) = 1;

%start with some baseline thresholding
mask(log10(hCI_mag/max(hCI_mag(:)))>-0.6) = 1; %-0.6 for Dragon Skin, -1.4 for ligs
%mask(log10(hires_mag/max(hires_mag(:)))>-0.6) = 1; %-0.6 for Dragon Skin, -1.4 for ligs
if plView3Ds, view3dgui(mask); end
if largestConnectedObjectMasking
    %view3dgui(mask);
    r = regionprops(mask>0);
    bwc = bwconncomp(mask>0);
    %find the largest connected object
    for k=1:length(r)
        area(k) = r(k).Area;
    end
    [MA, idx] = max(area);
    mask = mask - mask; %makes it zeros again;
    %only include that object
    mask(bwc.PixelIdxList{idx}) = 1; mask(mask==0) = nan;
end


%additional masking of other regions, for example
%mask([1:79 155:end],:,:) = nan;
if plView3Ds
    view3dgui(mask);
end
%mask(hCI_mag>0) = 1;

hCIx = interp3(cx,cy,cz,data(:,:,:,1),hx,hy,hz);
hCIy = interp3(cx,cy,cz,data(:,:,:,2),hx,hy,hz);
hCIz = interp3(cx,cy,cz,data(:,:,:,3),hx,hy,hz);

figure(31);
subplot(1,2,1);
imagesc(nanmean(real(hCIx),3)); colormap cividis; axis image;
set(gca,'Ydir','normal');
subplot(1,2,2);
imagesc(nanmean(imag(hCIx),3)); colormap cividis; axis image;
set(gca,'Ydir','normal');

shp = alphaShape(hx(mask==1),hy(mask==1),hz(mask==1))
%figure; plot(shp);

%Fix the mask to remove holes using the alphashape
TF = inShape(shp,hx,hy,hz); mask_ = nan(size(TF)); mask_(TF) = 1;
mask = mask_;


%% Check if there's a mask

%try load([fnum '_ref_vol_mask.mat'])
if callanmask
    try load('vol_mask.mat')
        mask = nan(size(tf));
        mask(tf) = 1;
        shp = alphaShape(hx(mask==1),hy(mask==1),hz(mask==1))
    end
end

figure(1);
% hold off;
% histogram(log10(hires_mag/max(hires_mag(:))))
% hold on;
% histogram(log10(hCI_mag/max(hCI_mag(:))))
figure; histogram(log10(hCI_mag/max(hCI_mag(:))),-3:0.01:0)
hold on; histogram(log10(mask.*hCI_mag/max(hCI_mag(:))),-3:0.01:0);


%Temporarily removed this for strain shift calcs, 180810
%hCIx = hCIx.*mask; hCIy = hCIy.*mask; hCIz = hCIz.*mask;
HCI = cat(4,hCI_mag/max(hCI_mag(:)),hCI_mag/max(hCI_mag(:)),hires_mag/max(hires_mag(:)));

figure(111); subplot(1,2,1); scatter(real(hCIx(mask==1)), imag(hCIx(mask==1)),'.'); axis image; hold on;

if prefilt
    %     hCIx = ndnanfilter(real(hCIx),'hamming',[5 5 5])+ 1i*ndnanfilter(imag(hCIx),'hamming',[5 5 5]);
    %     hCIy = ndnanfilter(real(hCIy),'hamming',[5 5 5])+ 1i*ndnanfilter(imag(hCIy),'hamming',[5 5 5]);
    %     hCIz = ndnanfilter(real(hCIz),'hamming',[5 5 5])+ 1i*ndnanfilter(imag(hCIz),'hamming',[5 5 5]);
    hCIx = ndnanfilter(hCIx,'hamming',[3 6 6]);
    hCIy = ndnanfilter(hCIy,'hamming',[3 6 6]);
    hCIz = ndnanfilter(hCIz,'hamming',[3 6 6]);
    
end

subplot(1,2,2); scatter(real(hCIx(mask==1)), imag(hCIx(mask==1)),'.'); axis image;

complexDivideUnwrap;
figure;
imagesc(squeeze(((strainimage(:,:,z,1,1)))).*squeeze(mask(:,:,z))); set(gca,'ydir','normal');
axis image; colormap magma; truesize;
caxis([-0.2 0.2]);

if saving
    save(['CI_', fnum, tag, '.mat'], 'hCIx', 'hCIy', 'hCIz');
end

figure(10); image(squeeze(HCI(:,:,z,:))); axis image;

hCI_phx = angle(hCIx);
hCI_phy = angle(hCIy);
hCI_phz = angle(hCIz);
%figure; view3dgui(hCI_phx)
%view3dgui(hCI_phy)

szI = size(HCI);


%% Run Julia script here to unwrap the 3D phase information
%
% load([fnum '_UnwrapX' tag '.mat']);
% load([fnum '_UnwrapY' tag '.mat']);
% load([fnum '_UnwrapZ' tag '.mat']);

%% Use matlab version of displacement calculation

if weighting
    w = ones(size(hCI_phx));
    w = w.*(hCI_mag/max(hCI_mag(:))).^2;
else
    w = ones(size(hCI_phx));
end
%try weighting by the phase information alone first

tic
juhCIx = WeightedUnwrap3(hCI_phx.*mask,w,mask);
toc
tic
juhCIy = WeightedUnwrap3(hCI_phy.*mask,w,mask);
toc
tic
juhCIz = WeightedUnwrap3(hCI_phz.*mask,w,mask);
toc

if saving
    save([fnum, '_UnwrapX' tag, '.mat'],'juhCIx');
    save([fnum, '_UnwrapY' tag, '.mat'],'juhCIy');
    save([fnum, '_UnwrapZ' tag, '.mat'],'juhCIz');
end

%%
%Displacement components u_i
ui{1} = juhCIx*lambda(1)/(2*pi);
ui{2} = juhCIy*lambda(2)/(2*pi);
ui{3} = juhCIz*lambda(3)/(2*pi);




figure(130);
subplot(2,3,1);
imagesc(-ui{1}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
title('u_1');
subplot(2,3,2);
imagesc(-ui{2}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
title('u_2');
subplot(2,3,3);
imagesc(-ui{3}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
title('u_3');


numfilts = 1;

for i=3:-1:1
    tu{i} = ui{i}(2:(end-1),2:(end-1),2:(end-1)); tu{i} = padarray(tu{i},[1 1 1],'symmetric');
    fu{i} = tu{i};
end

if ufilt
    for k = 1:numfilts
        %fu1 = imgaussfilt3(fu1,1,'FilterDomain','frequency');
        %fu2 = imgaussfilt3(fu2,1,'FilterDomain','frequency');
        %fu3 = imgaussfilt3(fu3,1,'FilterDomain','frequency');
        for i=3:-1:1
            %fu{i} = imgaussfilt3(fu{i});
            %fu{i} = imguidedfilter(fu{i});
            fu{i} = ndnanfilter(ui{i},'rectwin',[3 3 3]);
            %fu{i} = ndnanfilter(ui{i},'hamming',[5 5 5]);
        end
    end
end

subplot(2,3,4);
imagesc(-fu{1}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
title('u_1');
subplot(2,3,5);
imagesc(-fu{2}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
title('u_2');
subplot(2,3,6);
imagesc(-fu{3}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
title('u_3');


%% 3D figures


if pl3d
    tic
    comp = 1;
    %For nicked:
    %using [-0.05 2.4] for u1 (though -1.2 for each component on unfilt), [-0.132 0.055] for u2, [-0.25 0.05] for u3;
    %For dogbone
    %filt: using [-2.25 0] for u1, [-0.070 0.035] for u2, [-0.070 0.035] for u3
    uLims =  [min(-ui{comp}(:).*mask(:))-lambda(comp)/10,max(-ui{comp}(:).*mask(:))+lambda(comp)/10];
    figure(300+comp);
    hold on; sh = plot(shp); sh.FaceAlpha = 0.2; sh.EdgeAlpha = 0.02; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
    hold on;
    numpatches = 21;
    %pvals = linspace(0,wiggle(1),numpatches);
    %pvals = linspace(min(-ui{1}(:).*mask(:))*1.1,max(-ui{1}(:).*mask(:))*1.1,numpatches);
    pvals = uLims(1):lambda(comp)/20:uLims(2);%linspace(uLims(1),uLims(2),numpatches);
    magmaVals = colormap(magma);
    uj = ui{comp}.*mask;
    for np = 1:length(pvals)%numpatches
        p = patch(isosurface(hx,hy,hz,-uj,pvals(np)));
        %p.FaceColor = magmaVals(round(np/numpatches*length(magmaVals)),:);
        p.FaceColor = magmaVals(round(np/length(pvals)*length(magmaVals)),:);
        p.EdgeColor = 'none';
        p.FaceAlpha = 0.4;
    end
    daspect([1 1 1]);
    campos(cp);
    xlim([-5 5]); ylim([-30 30]); zlim([-5 5]);
    hold off;
    toc
    
    if ufilt
        tic
        figure(400+comp);
        numpatches = 21;
        %pvals = linspace(min(-fu{1}(:).*mask(:))*1.1,max(-fu{1}(:).*mask(:))*1.1,numpatches);
        %pvals = linspace(0,wiggle(1),numpatches);
        hold on; sh = plot(shp); sh.FaceAlpha = 0.2; sh.EdgeAlpha = 0.02; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
        pvals = linspace(-uLims(2),-uLims(1),numpatches);
        magmaVals = colormap(magma);
        fu1 = fu{1}.*mask;
        for np = 1:numpatches
            p = patch(isosurface(hx,hy,hz,-fu1,pvals(np)));
            p.FaceColor = magmaVals(round(np/numpatches*length(magmaVals)),:);
            p.EdgeColor = 'none';
            p.FaceAlpha = 0.4;
        end
        daspect([1 1 1]);
        campos(cp);
        hold off;
        xlim([-5 5]); ylim([-30 20]); zlim([-5 5]);
        
        toc
    end
end

%make ui,j non-scaled, with gradientN running perhaps optimal-n tap filter
[u12n, u11n, u13n] = gradientN(-(fu{1}.*mask),'optimal5');
[u22n, u21n, u23n] = gradientN(-(fu{2}.*mask),'optimal5');
[u32n, u31n, u33n] = gradientN(-(fu{3}.*mask),'optimal5');


%for pick z:
plsc = [hsc(1) hsc(2)];

%Grad_j u_i == u_i,j
Gjui{1,2} = u12n/hsc(2); Gjui{1,1} = u11n/hsc(1); Gjui{1,3} = u13n/hsc(3);
Gjui{2,2} = u22n/hsc(2); Gjui{2,1} = u21n/hsc(1); Gjui{2,3} = u23n/hsc(3);
Gjui{3,2} = u32n/hsc(2); Gjui{3,1} = u31n/hsc(1); Gjui{3,3} = u33n/hsc(3);

%Lagrangian strain tensor E_ij
Eij = calculateEij(Gjui);

%Small strain tensor eps_ij
epsij = calculateEpsij(Gjui);

%Von mises strain E_m
Em = calculateEm(Eij);


if pl3d
    comp = 1;
    XYZC{1} = hx; XYZC{2} = hy; XYZC{3} = hz; XYZC{4} = -ui{comp}.*mask;% - (max(-ui{comp}(:).*mask(:))-min(-ui{comp}(:).*mask(:)))/2;
    f = scatterColor3(XYZC, 21, magmaVals, 0, 0.1, 1, cp,uLims);
    hold off;
    
    if ufilt
        comp = 1;
        XYZC{1} = hx; XYZC{2} = hy; XYZC{3} = hz; XYZC{4} = -fu{comp}.*mask;% - (max(-ui{comp}(:).*mask(:))-min(-ui{comp}(:).*mask(:)))/2;
        f = scatterColor3(XYZC, 21, magmaVals, 0, 0.1, 1, cp ,uLims);
        hold off;
    end
    
    %ij(1) = 2; ij(2) = 2;
    XYZC{1} = hx; XYZC{2} = hy; XYZC{3} = hz; XYZC{4} = Em.*mask;%Eij{ij(1),ij(2)}.*mask;% - (max(-ui{comp}(:).*mask(:))-min(-ui{comp}(:).*mask(:)))/2;
    f = scatterColor3(XYZC, 21, magmaVals, 1.4, 0.2, 1, cp, [0 0.2]);
    hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.02; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
    hold off;
    
    
    %     Egood = ~isnan(Em);
    %     ids = find(Egood);
    %     colors = ceil(Em/max(Em(:))*64);
    %     C1 = nan(size(colors)); C2 = C1; C3 = C1;
    %     for cc = 1:sum(Egood(:))
    %         C1(ids(cc)) = magmaVals(colors(ids(cc)),1);
    %         C2(ids(cc)) = magmaVals(colors(ids(cc)),2);
    %         C3(ids(cc)) = magmaVals(colors(ids(cc)),3);
    %     end
    %     figure(202);
    %     hold on;
    %     numtiers = 20;
    %     gamma = 1.4;
    %     maxalpha = 0.8;
    %     for nt = 1:numtiers
    %     idxs = find(Egood & colors>((nt-1)/numtiers*64) & colors<=(nt/numtiers*64));
    %     %f2{nt} = scatter3(hx(Egood),hy(Egood),hz(Egood),1,cat(2,C1(Egood),C2(Egood),C3(Egood)));
    %     f2{nt} = scatter3(hx(idxs),hy(idxs),hz(idxs),2,cat(2,C1(idxs),C2(idxs),C3(idxs)));
    %     f2{nt}.MarkerEdgeAlpha = (nt/numtiers)^gamma*maxalpha;
    %     f2{nt}.MarkerFaceAlpha = (nt/numtiers)^gamma*maxalpha;
    %     end
    %     daspect([1 1 1]);
    %     campos([90 180 105]);
    
end

%%

%Lagrangian strain figure
%Normal components
figure(131);
subplot(2,3,1);
imagesc(Eij{1,1}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
caxis([-0.1 0.1]);
title('E_{11}');
subplot(2,3,2);
imagesc(Eij{2,2}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
caxis([-0.1 0.1]);
title('E_{22}');
subplot(2,3,3);
imagesc(Eij{3,3}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
caxis([-0.1 0.1]);
title('E_{33}');
%Shear components
subplot(2,3,4);
imagesc(Eij{1,2}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
caxis([-0.1 0.1]);
title('E_{12}');
subplot(2,3,5);
imagesc(Eij{2,3}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
caxis([-0.1 0.1]);
title('E_{23}');
subplot(2,3,6);
imagesc(Eij{1,3}(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
caxis([-0.1 0.1]);
title('E_{13}');

figure(132);
imagesc(Em(:,:,z).*mask(:,:,z)); set(gca,'ydir','normal'); axis image; colormap magma;
title('E_{m}');


%
% figure(20); histogram(reshape(Eij{1,1}(125:175,:),[1 numel(Eij{1,1}(125:175,:))])); hold on;
% histogram(reshape(Eij{2,2}(125:175,:),[1 numel(Eij{2,2}(125:175,:))]));
% histogram(reshape(Eij{3,3}(125:175,:),[1 numel(Eij{3,3}(125:175,:))]));
%
% gageE11 = Eij{1,1}(125:175,:,:).*mask(125:175,:,:);
% nonnan = sum(~isnan(gageE11(:)));
% gageE22 = Eij{2,2}(125:175,:,:).*mask(125:175,:,:);
% gageE33 = Eij{3,3}(125:175,:,:).*mask(125:175,:,:);
%
%
%
% nanmedian(gageE11(:))
% nanstd(gageE11(:))
%
% nanmedian(gageE22(:))
% nanstd(gageE22(:))
%
% nanmedian(gageE33(:))
% nanstd(gageE33(:))


figure(10);
subplot(1,5,1);
image(squeeze(HCI(:,:,z,:))); axis image; set(gca,'ydir','normal'); colorbar;
subplot(1,5,2);
imagesc(mask(:,:,z)); axis image; set(gca,'ydir','normal'); colorbar;
subplot(1,5,3);
imagesc(squeeze(-juhCIx(:,:,z)).*mask(:,:,z)); axis image; set(gca,'ydir','normal'); colormap magma;
colorbar;
%caxis([-2*pi 10*pi])
%subplot(1,5,3);
% medu1 = nanmedian(ui{1}(:).*mask(:));
% imagesc(squeeze(-(ui{1}(:,:,z)).*mask(:,:,z))); axis image; set(gca,'ydir','normal'); colormap magma;
% colorbar;
% caxis([-0.5 2.5]);
subplot(1,5,4);
medfu1 = nanmedian(fu{1}(:).*mask(:));
imagesc(squeeze(-(fu{1}(:,:,z)).*mask(:,:,z))); axis image; set(gca,'ydir','normal'); colormap magma;
colorbar;
%caxis([-0.5 1.5]);
subplot(1,5,5);
%imAlpha=ones(size(squeeze(Eij{1,1}(:,:,z))));
%imAlpha(isnan(Eij{1,1}(:,:,z)))=0;
%imagesc(squeeze(Eij{1,1}(:,:,z).*mask(:,:,z)),'AlphaData',imAlpha); axis image; set(gca,'ydir','normal');
imagesc(squeeze(Eij{1,1}(:,:,z).*mask(:,:,z))); axis image; set(gca,'ydir','normal');
%colormap(flipud([brewermap(256,'RdBu'); 0 0 0]));
colorbar;
%caxis([-0.1 0.1]);

if saving
    save([fnum '_mechvars.mat'], 'ui','Eij','epsij','Gjui','mask');
end

%figure(130); umax = max(sqrt(ux(:).^2+uy(:).^2+uz(:).^2));
%image(cat(3,ux/max(ux(:)),uy/max(uy(:)),uz/max(uz(:))))

