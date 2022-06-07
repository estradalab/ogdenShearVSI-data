%This script is the run function for the image processing related to the
%ligament stretching project.

%Written by Jon Estrada

clear; %close all;
%Uncomment to get rid of the view3dgui plots

% Use this to create DESTE file
% e.g. og=ogden_deste_3d(1958,1904,1931,1836,'blurvec',[0.8,0.8,0.8]);
%                         PE   SL   RO   Rf
% For 20211012-ogdenss

try
    delete(findall(0));
catch
end

sampleName = '20211012-ogdenss';
% sampleName = '20220124-ogdenss_5MMH_apod_64_16'; % Most optimal
% sampleName = '20220209-ogdenss_5MMH_translation';
% sampleName = '20220209-ogdenss_2moreloadsteps';
% sampleName = '20220209-ogdenss_translation';  
cd(sampleName)

switch sampleName
    case '20211012-ogdenss'
        fnums = {'1931'};
        refnum = '1836';
        fixedpt = [90,5,5];
        maskparam = 750;
    case '20220124-ogdenss_5MMH'
        % fnums = {'0924'}; % 1 Load step (5 mm)
        fnums = {'1139'}; % 2.5mm
        refnum = '0852';
        fixedpt = [101,9,15];
    case '20220124-ogdenss_5MMH_diffref'
        % fnums = {'1139','0924'}; % 2 Load steps
        fnums = {'0924'}; % 1 Load step (5 mm)
        refnum = '0843';
        fixedpt = [101,9,15];
    case '20220124-ogdenss_5MMH_2.5mm_lam2'
        fnums = {'1139'}; % 2.5 mm displacement
        refnum = '0852';
        fixedpt = [101,9,15];
    case '20220124-ogdenss_5MMH_2loadsteps'
        fnums = {'1139','0924'};
        refnum = '0852';
        fixedpt = [101,9,15];
    case '20220124-ogdenss_5MMH_apod_64_16' % Most optimal
        fnums = {'1139','0924'}; % 2 load steps
        refnum = '0852';
        fixedpt = [101,9,15];
        maskparam = 1300;
    case '20220124-ogdenss_5MMH_apod_128_32'
        fnums = {'0924'}; % 5mm
        refnum = '0852';
        fixedpt = [101,9,15];
    case '20220124-ogdenss_5MMH_apod_8_2'
        fnums = {'0924'}; % 5mm
        refnum = '0852';
        fixedpt = [101,9,15];
    case '20220209-ogdenss_5MMH_translation'
        fnums = {'1114'}; % 5mm translation
        refnum = '1215';
        fixedpt = [101,9,15];
        maskparam = 850;
    case '20220209-ogdenss_2moreloadsteps'
        fnums = {'1456','1336'};
        refnum = '1554';
        fixedpt = [89,6,6];
        maskparam = 600;
    case '20220209-ogdenss_translation'
        fnums = {'1657'};
        refnum = '1731';
        fixedpt = [89,6,6];
        maskparam = 600;
end

% hfilts = {[0 0 0],[1 1 1],[2 2 2],[3 3 3],[4 4 4],[8 8 8]};
% hfilts = {[0 0 0],[1 1 1],[2 2 2]};
hfilts = {[2 2 2]};
for hidx = 1:length(hfilts)
    for idx = 1:length(fnums)
        %         try
        %             delete(findall(0));
        %         catch
        %         end
        
        %%%% User-Defined Inputs %%%%
        %number of the reference image as a string
        fnum = fnums{idx};
        %tag for output files, e.g. 'nofilt', 'prekfilt', 'postfilt'
        tag = 'CD';
        %camera position/angle for the 3d plots
        %cp = [91.7068427025730,2.63079261003728,-142.224656664244];%
        cp = [90 180 105]; %cp = [40.1324   -1.2602 -145.6844];
        %chosen z value for plots
        z = 48; %48 for ligaments, 63 for dogbone
        
        %%%% Toggle on/off Options %%%%
        %save unwrap files, mechanical fields
        saving = false;
        %view 3d plots
        pl3d = true;
        %use origdata instead of HIRES
        origTF = true;
        %If we want to Gaussian filter (blur3d) complex data before unwrapping
        prefilt = false;
        hfilt = hfilts{hidx};
        postfilt = true;
        %Simple mask threshold on a log scale
        try
            maskthresh = -maskparam/1000;
        catch
            maskthresh = -0.6; % Generally for silicone
        end
        reps = 0;
        %filter u after unwrapping
        ufilt = false;
        %use weighted reliability algorithm, which fills in edges last
        weighting = true;
        %only include the largest single connected mask object
        largestConnectedObjectMasking=false;
        %open view3dgui plots
        plView3Ds = true;
        %Dilate mask
        dilatemask = false;
        fillmiss = true;
        %Corrected MRI
        corrMRI = true;
        
        
        %%%%%%%%
        
        hstr = [num2str(hfilt(1)) num2str(hfilt(2)) num2str(hfilt(3))];
        
        %loads the MR info
        fn = glob(['DESTE*' fnum '*.mat']);
        reffn = glob(['DESTE*' refnum '*.mat']);
        load(fn{1});
        L = lambda;
        try
            hsc = [abs(HIRES.axis1(2)-HIRES.axis1(1)), abs(HIRES.axis2(2)-HIRES.axis2(1)), abs(HIRES.axis3(2)-HIRES.axis3(1))];
            osc = [abs(origpars(1).axis1(2)-origpars(1).axis1(1)), abs(origpars(1).axis2(2)-origpars(1).axis2(1)), abs(origpars(1).axis3(2)-origpars(1).axis3(1))];
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
            hsc = [abs(HIRES.axis1(2)-HIRES.axis1(1)), abs(HIRES.axis2(2)-HIRES.axis2(1)), abs(HIRES.axis3(2)-HIRES.axis3(1))];
        end
        %reload to overwrite any vars with the same name
        load(fn{1});
        
        %Make a mask out of the HIRES info
        if ~origTF
            anat_mag = HIRES.magnitude;%flip(HIRES.magnitude,3);
            [ax,ay,az] = meshgrid(HIRES.axis2,HIRES.axis1,HIRES.axis3);
            try
                CI_mag = (abs(data(:,:,:,1))+abs(data(:,:,:,2))+abs(data(:,:,:,3)))/3;
            catch
                data = cat(4,data,data,data);
                CI_mag = (abs(data(:,:,:,1))+abs(data(:,:,:,2))+abs(data(:,:,:,3)))/3;
            end
        else
            origvol = load(reffn{1},'origdata','origpars');
            anat_mag = sqrt(max(abs(origvol.origdata(:,:,:,1)),4));
            [ax,ay,az] = meshgrid(origvol.origpars(1).axis2,origvol.origpars(1).axis1,origvol.origpars(1).axis3);
            CI_mag = (abs(origvol.origdata(:,:,:,1))+abs(origvol.origdata(:,:,:,2))+abs(origvol.origdata(:,:,:,3)))/3;
            data = origdata;
            axis1 = origpars(1).axis1; axis2 = origpars(1).axis2; axis3 = origpars(1).axis3;
        end
        [cx,cy,cz] = meshgrid(axis2,axis1,axis3);
        
        if corrMRI
            load('allAlphas.mat','final_a1fun','final_a2fun','final_a3fun',...
                'final_a11fun','final_a12fun','final_a13fun', ...
                'final_a21fun','final_a22fun','final_a23fun', ...
                'final_a31fun','final_a32fun','final_a33fun');
            corrGU{1,1} = final_a11fun(cy,cx,cz);
            corrGU{1,2} = final_a12fun(cy,cx,cz);
            corrGU{1,3} = final_a13fun(cy,cx,cz);
            corrGU{2,1} = final_a21fun(cy,cx,cz);
            corrGU{2,2} = final_a22fun(cy,cx,cz);
            corrGU{2,3} = final_a23fun(cy,cx,cz);
            corrGU{3,1} = final_a31fun(cy,cx,cz);
            corrGU{3,2} = final_a32fun(cy,cx,cz);
            corrGU{3,3} = final_a33fun(cy,cx,cz);
            
            corrU{1} = final_a1fun(cy,cx,cz);
            corrU{2} = final_a2fun(cy,cx,cz);
            corrU{3} = final_a3fun(cy,cx,cz);
            
        end
        
        
        hCI_mag = interp3(cx,cy,cz,CI_mag,ax,ay,az);

%% Thresholding
%         mask = nan(size(hCI_mag)); %mask(hires_mag>2E4) = 1;
%         %start with some baseline thresholding
%         mask(log10(hCI_mag/max(hCI_mag(:)))>maskthresh) = 1; %e.g. -0.6 for Dragon Skin, -1.4 for ligs\
%         range = [1:7 46:63 102:128];
%         mask(range,:,:)=NaN; % Plate for 20220124-ogdenss_5MMH
%         mask(:,:,31:32)=NaN; % Noise at end for 20220124-ogdenss_5MMH
        
%% Physical Sizes
        switch sampleName
            case '20211012-ogdenss'
                ranget = [1:16 55:71 110:128];  % 38/17/38 px thickness
                rangew = [1:3 30:32];           % 26 px width
                rangel = [1:4 31:32];           % 26 px length
                mask = ones(size(hCI_mag));
                mask(ranget,:,:) = NaN;
                mask(:,rangew,:) = NaN;
                mask(:,:,rangel) = NaN;
                % mask(log10(hCI_mag/max(hCI_mag(:)))<maskthresh) = NaN;
            case '20220124-ogdenss_5MMH_apod_64_16' % Most optimal
                ranget = [1:12 49:65 102:128];  % 36 px thickness
                rangew = [1:3 30:32];           % 27 px width
                rangel = [1:3 30:32];           % 27 px length
                mask = ones(size(hCI_mag));
                mask(ranget,:,:) = NaN;
                mask(:,rangew,:) = NaN;
                mask(:,:,rangel) = NaN;
                mask(log10(hCI_mag/max(hCI_mag(:)))<maskthresh) = NaN;
                for k = 1:4
                    if k == 1
                        mask(14,:,:) = mask(15,:,:);
                        mask(46,:,:) = mask(45,:,:);
                        mask(68,:,:) = mask(67,:,:);
                        mask(98,:,:) = mask(97,:,:);
                    elseif k == 2
                        mask(13,:,:) = mask(15,:,:);
                        mask(47,:,:) = mask(45,:,:);
                        mask(67,:,:) = mask(67,:,:);
                        mask(99,:,:) = mask(97,:,:);
                    elseif k == 3
                        mask(48,:,:) = mask(45,:,:);
                        mask(66,:,:) = mask(67,:,:);
                        mask(100,:,:) = mask(97,:,:);
                    elseif k == 4
                        mask(101,:,:) = mask(97,:,:);
                    end
                end
                mask([13:48 66:101],[4:11 22:29],[4:29]) = ones(size(mask([13:48 66:101],[4:11 22:29],[4:29])));
                mask([13:48 66:101],[4:29],[4:8 25:29]) = ones(size(mask([13:48 66:101],[4:29],[4:8 25:29])));

            case '20220209-ogdenss_5MMH_translation'
                ranget = [1:15 54:71 103:128];  % 35 px thickness
                rangew = [1:3 31:32];           % 26 px width
                rangel = [1:3 31:32];           % 26 px length
                mask = ones(size(hCI_mag));
                mask(ranget,:,:) = NaN;
                mask(:,rangew,:) = NaN;
                mask(:,:,rangel) = NaN;
                mask(log10(hCI_mag/max(hCI_mag(:)))<maskthresh) = NaN;
            case '20220209-ogdenss_2moreloadsteps'
                ranget = [1:15 54:70 109:128];  % 38/17/38 px thickness
                rangew = [1:4 31:32];           % 26 px width
                rangel = [1:4 31:32];           % 26 px length
                mask = ones(size(hCI_mag));
                mask(ranget,:,:) = NaN;
                mask(:,rangew,:) = NaN;
                mask(:,:,rangel) = NaN;
                % mask(log10(hCI_mag/max(hCI_mag(:)))<maskthresh) = NaN;
            case '20220209-ogdenss_translation'
                ranget = [1:13 50:66 103:128];  % 35 px thickness
                rangew = [1:4 31:32];           % 25 px width
                rangel = [1:3 30:32];           % 25 px length
                mask = ones(size(hCI_mag));
                mask(ranget,:,:) = NaN;
                mask(:,rangew,:) = NaN;
                mask(:,:,rangel) = NaN;
            otherwise
                mask = nan(size(hCI_mag)); %mask(hires_mag>2E4) = 1;
                %start with some baseline thresholding
                mask(log10(hCI_mag/max(hCI_mag(:)))>maskthresh) = 1; %e.g. -0.6 for Dragon Skin, -1.4 for ligs\
        end

        %mask(log10(hires_mag/max(hires_mag(:)))>-0.6) = 1; %-0.6 for Dragon Skin, -1.4 for ligs
        if plView3Ds
            %view3dgui(mask); 
        end
        if largestConnectedObjectMasking
            %view3dgui(mask);
            r = regionprops(mask>0);
            bwc = bwconncomp(mask>0);
            %find the largest connected object...
            for k=1:length(r)
                area(k) = r(k).Area;
            end
            [MA, idxz] = max(area);
            mask = mask - mask; %makes it zeros again;
            %then remove all the non-connected smaller objects
            mask(bwc.PixelIdxList{idxz}) = 1; mask(mask==0) = nan;
            %         try
            %             %    load('mask.mat');
            %             %catch
            %             load('mask_HIRES.mat');
            %             mask = round(mask_HIRES(:,1:2:127,:)/2+mask_HIRES(:,2:2:128,:)/2);
            %         end
        end
        
        %Look at the mask if plot view 3D is on
        if plView3Ds
            % view3dgui(mask);
        end
        %mask(hCI_mag>0) = 1;
        
        if prefilt
            %Perform blurring filter based on the size of hfilt
%             datax = blur3d(data(:,:,:,1),'vox',hfilt);
%             datay = blur3d(data(:,:,:,2),'vox',hfilt);
%             dataz = blur3d(data(:,:,:,3),'vox',hfilt);

            datax = ndnanfilter(data(:,:,:,1),'hamming',hfilt);
            datay = ndnanfilter(data(:,:,:,2),'hamming',hfilt);
            dataz = ndnanfilter(data(:,:,:,3),'hamming',hfilt);
            
            % ndnanfilter(Fij{ii,jj},'hamming',hfilt).*mask;

            %If preferred, do a sequence of small filterings
            %Note that N filts is equivalent to a single filt of size sqrt(N)
            if reps ~=0
                for m=1:reps
                    datax = blur3d(datax,'vox',hfilt);
                    datay = blur3d(datay,'vox',hfilt);
                    dataz = blur3d(dataz,'vox',hfilt);
                end
            end
        else
            datax = data(:,:,:,1);
            datay = data(:,:,:,2);
            dataz = data(:,:,:,3);
        end
        
        hCIx = interp3(cx,cy,cz,datax,ax,ay,az);
        hCIy = interp3(cx,cy,cz,datay,ax,ay,az);
        hCIz = interp3(cx,cy,cz,dataz,ax,ay,az);
        
        %axmx = ndgrid(axis1,axis2,axis3);
        
        chy = 17;
%         figure(31);
%         subplot(1,4,1);
%         imagesc(axis3, axis1, squeeze(real(hCIx(:,chy,:)))); colormap cividis; axis image;
%         set(gca,'Ydir','normal');
%         subplot(1,4,2);
%         imagesc(axis3, axis1, squeeze(imag(hCIx(:,chy,:)))); colormap cividis; axis image;
%         set(gca,'Ydir','normal');
%         subplot(1,4,3);
%         imagesc(axis3, axis1, squeeze(abs(hCIx(:,chy,:)))); colormap cividis; axis image;
%         set(gca,'Ydir','normal');
%         subplot(1,4,4);
%         imagesc(axis3, axis1, squeeze(angle(hCIx(:,chy,:).*mask(:,chy,:)))); colormap cividis; axis image;
%         set(gca,'Ydir','normal');
        
        %Adding in a radius, 211030
        shp = alphaShape(ax(mask==1),ay(mask==1),az(mask==1),1);
        %figure; plot(shp);
        
        %Fix the mask to remove internal holes using the alphashape
        TF = inShape(shp,ax,ay,az); mask_ = nan(size(TF)); mask_(TF) = 1;
        mask = mask_;
        % view3dgui(mask_)
        % view3dgui(mask)
        
        %% Check if there's a separate mask
        
        %  try load('vol_mask.mat')%load([fnum '_ref_vol_mask.mat'])
        %      mask = nan(size(tf));
        %      mask(tf) = 1;
        %      shp = alphaShape(hx(mask==1),hy(mask==1),hz(mask==1))
        %  end
        
        %Signal histogram from the MR - blue = all points, orange = those in mask
        figure; histogram(log10(hCI_mag/max(hCI_mag(:))),-3:0.01:0);
        hold on; histogram(log10(mask.*hCI_mag/max(hCI_mag(:))),-3:0.01:0);
        
        
        %Temporarily removed this for strain shift calcs, 180810
        %hCIx = hCIx.*mask; hCIy = hCIy.*mask; hCIz = hCIz.*mask;
        HCI = cat(4,hCI_mag/max(hCI_mag(:)),hCI_mag/max(hCI_mag(:)),anat_mag/max(anat_mag(:)));
        
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
        
        LResMask = interp3(ax,ay,az,mask,cx,cy,cz);
        % filts{1} = 'optimal3'; %filts{2} = 'optimal3'; %filts{3} = 'optimal5';
        filts{1} = 'fb';
        % filts{1} = 'sobel';
        
        mask_ = mask>0;
        dmask = imdilate(mask_, ones(5,5,5));
        nmask = nan(size(mask)); nmask(dmask) = 1;
        
        % minimask = mask(101:150,:,:);
        % XApx = sum(~isnan(minimask(:))/50);
        % XA = XApx*hsc(2)*hsc(3)/1E6;
        
        for j = 1:length(filts)
            
            filt = filts{j};
%             if dilatemask
%                 [Qij{1,1}(:,:,:), Qij{1,2}(:,:,:), Qij{1,3}(:,:,:)] = complexDivideN(hCIx.*nmask,filt);
%                 [Qij{2,1}(:,:,:), Qij{2,2}(:,:,:), Qij{2,3}(:,:,:)] = complexDivideN(hCIy.*nmask,filt);
%                 [Qij{3,1}(:,:,:), Qij{3,2}(:,:,:), Qij{3,3}(:,:,:)] = complexDivideN(hCIz.*nmask,filt);
%             else
%                 [Qij{1,1}(:,:,:), Qij{1,2}(:,:,:), Qij{1,3}(:,:,:)] = complexDivideN(hCIx.*mask,filt);
%                 [Qij{2,1}(:,:,:), Qij{2,2}(:,:,:), Qij{2,3}(:,:,:)] = complexDivideN(hCIy.*mask,filt);
%                 [Qij{3,1}(:,:,:), Qij{3,2}(:,:,:), Qij{3,3}(:,:,:)] = complexDivideN(hCIz.*mask,filt);
%             end

            if dilatemask
                %Qij{1,1} = padarray(angle(divv(hCIx.*nmask./abs(hCIx),1)/(origpars(1).axis1(2)-origpars(1).axis1(1))),[1 0 0],'replicate','both'); 
                %Qij{1,2} = padarray(angle(divv(hCIx.*nmask./abs(hCIx),2)/(origpars(1).axis2(2)-origpars(1).axis2(1))),[0 1 0],'replicate','both');
                %Qij{1,3} = padarray(angle(divv(hCIx.*nmask./abs(hCIx),3)/(origpars(1).axis3(2)-origpars(1).axis3(1))),[0 0 1],'replicate','both');
                %Qij{2,1} = padarray(angle(divv(hCIy.*nmask./abs(hCIy),1)/(origpars(1).axis1(2)-origpars(1).axis1(1))),[1 0 0],'replicate','both');
                %Qij{2,2} = padarray(angle(divv(hCIy.*nmask./abs(hCIy),2)/(origpars(1).axis2(2)-origpars(1).axis2(1))),[0 1 0],'replicate','both');
                %Qij{2,3} = padarray(angle(divv(hCIy.*nmask./abs(hCIy),3)/(origpars(1).axis3(2)-origpars(1).axis3(1))),[0 0 1],'replicate','both');
                %Qij{3,1} = padarray(angle(divv(hCIz.*nmask./abs(hCIz),1)/(origpars(1).axis1(2)-origpars(1).axis1(1))),[1 0 0],'replicate','both');
                %Qij{3,2} = padarray(angle(divv(hCIz.*nmask./abs(hCIz),2)/(origpars(1).axis2(2)-origpars(1).axis2(1))),[0 1 0],'replicate','both');
                %Qij{3,3} = padarray(angle(divv(hCIz.*nmask./abs(hCIz),3)/(origpars(1).axis3(2)-origpars(1).axis3(1))),[0 0 1],'replicate','both');
                Qij{1,1} = padarray(angle(divv(hCIx.*nmask./abs(hCIx),1)),[1 0 0],'replicate','both'); 
                Qij{1,2} = padarray(angle(divv(hCIx.*nmask./abs(hCIx),2)),[0 1 0],'replicate','both');
                Qij{1,3} = padarray(angle(divv(hCIx.*nmask./abs(hCIx),3)),[0 0 1],'replicate','both');
                Qij{2,1} = padarray(angle(divv(hCIy.*nmask./abs(hCIy),1)),[1 0 0],'replicate','both');
                Qij{2,2} = padarray(angle(divv(hCIy.*nmask./abs(hCIy),2)),[0 1 0],'replicate','both');
                Qij{2,3} = padarray(angle(divv(hCIy.*nmask./abs(hCIy),3)),[0 0 1],'replicate','both');
                Qij{3,1} = padarray(angle(divv(hCIz.*nmask./abs(hCIz),1)),[1 0 0],'replicate','both');
                Qij{3,2} = padarray(angle(divv(hCIz.*nmask./abs(hCIz),2)),[0 1 0],'replicate','both');
                Qij{3,3} = padarray(angle(divv(hCIz.*nmask./abs(hCIz),3)),[0 0 1],'replicate','both');
            else
%                 Qij{1,1} = padarray(angle(divv(hCIx.*mask./abs(hCIx),1)/(origpars(1).axis1(2)-origpars(1).axis1(1))),[1 0 0],'replicate','both');
%                 Qij{1,2} = padarray(angle(divv(hCIx.*mask./abs(hCIx),2)/(origpars(1).axis2(2)-origpars(1).axis2(1))),[0 1 0],'replicate','both');
%                 Qij{1,3} = padarray(angle(divv(hCIx.*mask./abs(hCIx),3)/(origpars(1).axis3(2)-origpars(1).axis3(1))),[0 0 1],'replicate','both');
%                 Qij{2,1} = padarray(angle(divv(hCIy.*mask./abs(hCIy),1)/(origpars(1).axis1(2)-origpars(1).axis1(1))),[1 0 0],'replicate','both');
%                 Qij{2,2} = padarray(angle(divv(hCIy.*mask./abs(hCIy),2)/(origpars(1).axis2(2)-origpars(1).axis2(1))),[0 1 0],'replicate','both');
%                 Qij{2,3} = padarray(angle(divv(hCIy.*mask./abs(hCIy),3)/(origpars(1).axis3(2)-origpars(1).axis3(1))),[0 0 1],'replicate','both');
%                 Qij{3,1} = padarray(angle(divv(hCIz.*mask./abs(hCIz),1)/(origpars(1).axis1(2)-origpars(1).axis1(1))),[1 0 0],'replicate','both');
%                 Qij{3,2} = padarray(angle(divv(hCIz.*mask./abs(hCIz),2)/(origpars(1).axis2(2)-origpars(1).axis2(1))),[0 1 0],'replicate','both');
%                 Qij{3,3} = padarray(angle(divv(hCIz.*mask./abs(hCIz),3)/(origpars(1).axis3(2)-origpars(1).axis3(1))),[0 0 1],'replicate','both');
                Qij{1,1} = padarray(angle(divv(hCIx.*mask./abs(hCIx),1)),[1 0 0],'replicate','both');
                Qij{1,2} = padarray(angle(divv(hCIx.*mask./abs(hCIx),2)),[0 1 0],'replicate','both');
                Qij{1,3} = padarray(angle(divv(hCIx.*mask./abs(hCIx),3)),[0 0 1],'replicate','both');
                Qij{2,1} = padarray(angle(divv(hCIy.*mask./abs(hCIy),1)),[1 0 0],'replicate','both');
                Qij{2,2} = padarray(angle(divv(hCIy.*mask./abs(hCIy),2)),[0 1 0],'replicate','both');
                Qij{2,3} = padarray(angle(divv(hCIy.*mask./abs(hCIy),3)),[0 0 1],'replicate','both');
                Qij{3,1} = padarray(angle(divv(hCIz.*mask./abs(hCIz),1)),[1 0 0],'replicate','both');
                Qij{3,2} = padarray(angle(divv(hCIz.*mask./abs(hCIz),2)),[0 1 0],'replicate','both');
                Qij{3,3} = padarray(angle(divv(hCIz.*mask./abs(hCIz),3)),[0 0 1],'replicate','both');
            end

            for i = 1:size(Qij{1,1},1)-1
                Qij_{1,1}(i,:,:) = (Qij{1,1}(i,:,:)+Qij{1,1}(i+1,:,:))/2;
                Qij_{2,1}(i,:,:) = (Qij{2,1}(i,:,:)+Qij{2,1}(i+1,:,:))/2;
                Qij_{3,1}(i,:,:) = (Qij{3,1}(i,:,:)+Qij{3,1}(i+1,:,:))/2;
            end
            for i = 1:size(Qij{2,2},2)-1
                Qij_{1,2}(:,i,:) = (Qij{1,2}(:,i,:)+Qij{1,2}(:,i+1,:))/2;
                Qij_{2,2}(:,i,:) = (Qij{2,2}(:,i,:)+Qij{2,2}(:,i+1,:))/2;
                Qij_{3,2}(:,i,:) = (Qij{3,2}(:,i,:)+Qij{3,2}(:,i+1,:))/2;
            end
            for i = 1:size(Qij{3,3},3)-1
                Qij_{1,3}(:,:,i) = (Qij{1,3}(:,:,i)+Qij{1,3}(:,:,i+1))/2;
                Qij_{2,3}(:,:,i) = (Qij{2,3}(:,:,i)+Qij{2,3}(:,:,i+1))/2;
                Qij_{3,3}(:,:,i) = (Qij{3,3}(:,:,i)+Qij{3,3}(:,:,i+1))/2;
            end
            clear Qij
            Qij = Qij_;
            
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
                        % Equation for 3 x 3 filter
                        % GUij{e,g}(:,:,:) = -Qij{e,g}(:,:,:)/(2*pi)*lambda(e)/(2*delta(g));
                        % Equation for 2 x 2 filter
                        GUij{e,g}(:,:,:) = -Qij{e,g}(:,:,:)/(2*pi)*lambda(e)/(delta(g));
                    catch
                        % Equation for 3 x 3 filter
                        % GUij{e,g}(:,:,:) = -Qij{e,g}(:,:,:)/(2*pi)*lambda/(2*delta(g));
                        % Equation for 2 x 2 filter
                        GUij{e,g}(:,:,:) = -Qij{e,g}(:,:,:)/(2*pi)*lambda/(delta(g));
                    end
                end
            end
            
            %     E11_err = -datax_err/(2*pi)*lambda(1)/(2*delta(1));
            %     E12_err = -datax_err/(2*pi)*lambda(1)/(2*delta(2));
            %     E22_err = -datay_err/(2*pi)*lambda(2)/(2*delta(2));
            %     E23_err = -datax_err/(2*pi)*lambda(2)/(2*delta(3));
            %     E13_err = -datax_err/(2*pi)*lambda(1)/(2*delta(3));
            %     E33_err = -dataz_err/(2*pi)*lambda(3)/(2*delta(3));
            dij = eye(3,3);
            for iii=1:3
                for jjj=1:3
                    Fij{iii,jjj} = GUij{iii,jjj}+dij(iii,jjj);
                end
            end
            
            tic
            Isz = size(Fij{1,1});
            F = cell(Isz);
            for ii=1:Isz(1)
                for jj=1:Isz(2)
                    for kk=1:Isz(3)
                        F_ = [Fij{1,1}(ii,jj,kk) Fij{1,2}(ii,jj,kk) Fij{1,3}(ii,jj,kk);...
                            Fij{2,1}(ii,jj,kk) Fij{2,2}(ii,jj,kk) Fij{2,3}(ii,jj,kk);...
                            Fij{3,1}(ii,jj,kk) Fij{3,2}(ii,jj,kk) Fij{3,3}(ii,jj,kk);];
                        F{ii,jj,kk} = F_;
                    end
                end
            end
            toc
            
            tic
            if fillmiss
                for ii=1:3
                    for jj=1:3
                        Fij_{ii,jj} = Fij{ii,jj};
                        %if ii==jj
                        %Try to suppress x-related duplication of noise
                        %by never making the 1-pad first
                        
                        %but also, always do the actual derivative
                        %direction pad last
                        %                             switch jj
                        %                                 case 1
                        %                                     order = [2 3 1];
                        %                                 case 2
                        %                                     order = [3 1 2];
                        %                                 case 3
                        %                                     order = [2 1 3];
                        %                             end
                        switch jj
                            case 1
                                order = [2 3 1];
                            case 2
                                order = [3 2 1];
                            case 3
                                order = [2 3 1];
                        end
                        %else
                        %    order = [6-(1+jj) jj 1];
                        %end
                        Fij{ii,jj} = fillmissing(Fij{ii,jj},'nearest',order(1));
                        Fij{ii,jj} = fillmissing(Fij{ii,jj},'nearest',order(2));
                        Fij{ii,jj} = fillmissing(Fij{ii,jj},'nearest',order(3));
                        %Changed the order to reverse - for 11 strain, do
                        %[3 2 1], for 23 strain do [1 2 3]
                        %Original code:
                        %order = [mod(jj,3) mod(jj+1,3) mod(jj+2,3)]; order(order==0) = 3;
                        %Fij{ii,jj} = fillmissing(Fij{ii,jj},'nearest',order(1));
                        %Fij{ii,jj} = fillmissing(Fij{ii,jj},'nearest',order(2));
                        %Fij{ii,jj} = fillmissing(Fij{ii,jj},'nearest',order(3));
                        sum(~isnan(Fij{1,1}(:)));
                        Fij{ii,jj} = Fij{ii,jj}.*mask;
                    end
                end
            end
            toc
            
%             if postfilt
%                 for ii=1:3
%                     for jj=1:3
%                         %Perform blurring filter based on the size of hfilt
%                         %Fij{ii,jj} = blur3d(Fij{ii,jj},'vox',hfilt);
%                         Fij{ii,jj} = ndnanfilter(Fij{ii,jj},'hamming',hfilt).*mask;
%                     end
%                 end
%             end
            
            for iii=1:3
                for jjj=1:3
                    GUij{iii,jjj} = Fij{iii,jjj} - dij(iii,jjj);
                end
            end
            
            %Calculate small strain eps, Lagrange-Green strain Eij, and effective strain Em
            epsij = calculateEpsij(GUij);
            Eij = calculateEij(GUij);
            Em = calculateEm(Eij);
            
            %             if ~prefilt
            %                 hstr = 'off';
            %             end
            
        end
        
        %Strain histogram
        figure(80); hold on; histogram(Eij{1,1});
        histogram(Eij{2,2}); histogram(Eij{3,3});
        
        %% Strain plots in 3D using my custom function scatterColor3
        
        if ~pl3d
            close all;
            %magmaVals = colormap(magma);
            rbVals = flipud(brewermap(256,'RdYlBu'));
            XYZC{1} = ax-median(ax(:)); XYZC{2} = ay-median(ay(:)); XYZC{3} = az-median(az(:));
            XYZC{4} = Eij{1,1}.*mask;
            %XYZC{4} = smooth3(Eij{1,1},'gaussian',[5 3 3]).*mask;%Eij{ij(1),ij(2)}.*mask;% - (max(-ui{comp}(:).*mask(:))-min(-ui{comp}(:).*mask(:)))/2;
            %XYZC{4} = ndnanfilter(Em,'hamming',[5 3 3]).*mask;
            
            falpha = 0.9;
            fxx = scatterColor3(XYZC, 21, rbVals, 1.0, falpha, 20, cp, [-0.25 0.25]);
            set(gca,'color','k');
            %hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.01; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
            %hold off;
            xlim([-5 5]); ylim([-25 25]); zlim([-5 5])
            %saveas(gcf,'wiggleE11.png')
            
            XYZC{4} = Eij{2,2}.*mask;
            fyy = scatterColor3(XYZC, 21, rbVals, 1.0, falpha, 20, cp, [-0.125 0.125]);
            set(gca,'color','k');
            %hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.01; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
            %hold off;
            xlim([-5 5]); ylim([-25 25]); zlim([-5 5])
            %saveas(gcf,'wiggleE22.png')
            
            XYZC{4} = Eij{3,3}.*mask;
            fzz = scatterColor3(XYZC, 21, rbVals, 1.0, falpha, 20, cp, [-0.125 0.125]);
            set(gca,'color','k');
            %hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.01; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
            %hold off;
            xlim([-5 5]); ylim([-25 25]); zlim([-5 5])
            %saveas(gcf,'wiggleE33.png')
            
            XYZC{4} = Eij{1,2}.*mask;
            fxy = scatterColor3(XYZC, 21, rbVals, 1.0, falpha, 20, cp, [-0.125 0.125]);
            set(gca,'color','k');
            xlim([-5 5]); ylim([-25 25]); zlim([-5 5])
            
            XYZC{4} = Eij{2,3}.*mask;
            fyz = scatterColor3(XYZC, 21, rbVals, 1.0, falpha, 20, cp, [-0.125 0.125]);
            xlim([-5 5]); ylim([-25 25]); zlim([-5 5])
            set(gca,'color','k');
            
            XYZC{4} = Eij{1,3}.*mask;
            fxz = scatterColor3(XYZC, 21, rbVals, 1.0, falpha, 20, cp, [-0.125 0.125]);
            set(gca,'color','k');
            xlim([-5 5]); ylim([-25 25]); zlim([-5 5])
            %hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.01; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
            %hold off;
            %xlim([-5 5]); ylim([-25 15]); zlim([0 10])
            %saveas(gcf,'wiggleE33.png')
            %
            %     rbVals = colormap(parula);
            %     XYZC{4} = Em.*mask;
            %     fm = scatterColor3(XYZC, 21, rbVals, 1.0, 0.6, 12, cp, [0 0.4]);
            %hold on; sh = plot(shp); sh.FaceAlpha = 0.6; sh.EdgeAlpha = 0.01; sh.FaceColor = 0.94*[1 1 1]; sh.LineWidth = 0.01;
            %hold off;
            %xlim([-5 5]); ylim([-25 15]); zlim([0 10])
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
        
        if postfilt && hfilt(1) ~= 0
            for ii=1:3
                for jj=1:3
                %Perform blurring filter based on the size of hfilt
                % Fij{ii,jj} = blur3d(Fij{ii,jj},'vox',hfilt);
                Fij{ii,jj} = ndnanfilter(Fij{ii,jj},'hamming',hfilt).*mask;
                end
                ui{ii} = ndnanfilter(ui{ii},'hamming',hfilt).*mask;
            end
        end
        
        if ~pl3d
            view3dgui(ui{1})
            view3dgui(ui{2})
            view3dgui(ui{3})
            
        end
        
        
        if corrMRI
            
            for iii=1:3
                for jjj=1:3
                    Fij{iii,jjj} = Fij{iii,jjj} - corrGU{iii,jjj}.*(ui{1}-ui{1}(fixedpt(1),fixedpt(2),fixedpt(3)));
                    GUij{iii,jjj} = Fij{iii,jjj} - dij(iii,jjj);
                end
                corrfact{iii} = (ui{1}-ui{1}(fixedpt(1),fixedpt(2),fixedpt(3))).*corrU{iii};
                corrfact{iii} = corrfact{iii}-corrfact{iii}(fixedpt(1),fixedpt(2),fixedpt(3));
                uicorr{iii} = ui{iii}-corrfact{iii};
                uicorr{iii} = 1/1.0634*(uicorr{iii}-uicorr{iii}(fixedpt(1),fixedpt(2),fixedpt(3))); %correcting for Uli's scale factor, 191204
            end
            
            %Calculate small strain eps, Lagrange-Green strain Eij, and effective strain Em
            epsij = calculateEpsij(GUij);
            Eij = calculateEij(GUij);
            Em = calculateEm(Eij);
            
            max(uicorr{1}(:))-min(uicorr{1}(:))
            max(ui{1}(:))-min(ui{1}(:))
        end
        
        F_t{idx} = Fij;
        for i = 1:3
            for j = 1:3
                E_t{idx}{i,j} = Eij{i,j}.*mask;
            end
        end
        prescribedU(idx) = wiggle(1);

        % Boundary condition
        uicorr{3}(1:64,:,:) = uicorr{3}(1:64,:,:) + (prescribedU(idx)-max(uicorr{3}(1:64,:,:),[],'all')); % x mm at traction surface
        uicorr{3}(64:128,:,:) = uicorr{3}(64:128,:,:) + (prescribedU(idx)-max(uicorr{3}(64:128,:,:),[],'all'));
        
        if corrMRI
            U_t{idx} = uicorr;
        else
            U_t{idx} = ui;
        end
        
        %Save into a file
        if saving
            save(['Fij_' filt 'h' hstr '_' fnum 'fm-newpad-ham-MRcorr-1last-220121.mat'], 'epsij', 'Eij', 'Fij', 'Em','mask', 'hfilt','dilatemask','fillmiss','uicorr','wiggle');
            %save(['MRI-3Ddefs_RectPrism_190919_' hstr 'hfilt_' fnum],'Fij','uicorr','mask','wiggle');
            
            %save(['MRI-3Ddefs_RectPrism_190919_' hstr 'hfilt'],'F_t','U_t','mask','prescribedU');
        end
        
        figure(32);
        subplot(1,3,1);
        imagesc(axis3, axis1, squeeze(U_t{idx}{1}(:,chy,:))); colormap cividis; axis image;
        set(gca,'Ydir','normal');
        subplot(1,3,2);
        imagesc(axis3, axis1, squeeze(U_t{idx}{2}(:,chy,:))); colormap cividis; axis image;
        set(gca,'Ydir','normal');
        title('Displacement field components, fixed Y-slice');
        subplot(1,3,3);
        imagesc(axis3, axis1, squeeze(U_t{idx}{3}(:,chy,:))); colormap cividis; axis image;
        set(gca,'Ydir','normal');
        
        
        
        chz = 13;
        figure(33);
        subplot(1,3,1);
        imagesc(axis2, axis1, squeeze(U_t{idx}{1}(:,:,chz))); colormap cividis; axis image;
        set(gca,'Ydir','normal');
        subplot(1,3,2);
        imagesc(axis2, axis1, squeeze(U_t{idx}{2}(:,:,chz))); colormap cividis; axis image;
        set(gca,'Ydir','normal');
        title('Displacement field components, fixed Z-slice');
        subplot(1,3,3);
        imagesc(axis2, axis1, squeeze(U_t{idx}{3}(:,:,chz))); colormap cividis; axis image;
        set(gca,'Ydir','normal');
        
        
        %calculatingAlpha
        %        imagesc(axis3, axis1, squeeze(angle(hCIx(:,chy,:).*mask(:,chy,:)))); colormap cividis; axis image;
        %        set(gca,'Ydir','normal');
        
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
        xlim([3 6]);
        ylim([3 6]);
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
    end
    
end
X = {ax;ay;az};

% Save mask_data_***.mat

save(['mask_data_' num2str(maskparam) '.mat'],'mask','mask_','maskthresh','osc');
save(['disp_data_' sampleName '.mat'],'F_t','E_t','U_t','mask','prescribedU','osc');
save(['refpositions_' sampleName '.mat'],'X');
cd ..