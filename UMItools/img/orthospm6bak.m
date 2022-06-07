function orthospm6(roi, threshold, onsets, window ,spm_file, anat_file, tseries_path, func_root, cx,cy,cz, voxelFile )
% function result = orthospm6(
%                         roi_radius,
%                         threshold , 
%                         onsets, 
%                         window 
%                         [,spm_file , 
%                          anat_file, 
%                          tseries_path, 
%                          func_root_name,
%			    x,y,z,
%                           voxelFile] )
%
%
% this function takes a file in analyze format and 
% displays orthogonal sections thru the planes interscting
% at x,y,z.  It then overlays a Statistical Parameter Map 
% thresholded at "threshold".
%
% the data extracted are those from the voxels specified in the 
% TEXT file "voxelFile".  (x,y,z) coordinates are just for display.
% The threshold doesn't matter.
%
% ** this function is NOT interactive on purpose **
% in fact, a lot of the input params are irrelevant...
%
% It saves the files:
% tdata.data - the time series extracted
% avg.dat - the averaged events
% pixels.dat - the xyz coords of the pixels used 
%  
% *** This version of the program calls event_avg.m 
% in order to average all the events into one.  No deconvolution!!
%
% onsets and window must be in scan units
 
global SPM_scale_factor 
%close all
figure
    if nargin < 5
        [name path] = uigetfile('*.img','Select SPM *.img file');
        spm_file = strcat(path,name)
    end
       
    % figure out name for spm file
    sz = size(spm_file);
    imgname = strcat(spm_file(1,1:sz(2)-4) , '.img');
    hdrname = strcat(spm_file(1,1:sz(2)-4) , '.hdr');
    %label=sprintf('FUNCTIONAL:  ....%s',imgname(end-15:end));
	
    spm_hdr = read_hdr(hdrname);
    spm_data = read_img2(spm_hdr,imgname);
    spm_data(find(spm_data==NaN))=0;
    spm_scale =  SPM_scale_factor;

    if nargin < 5
        [name path] = uigetfile('*.img','Select anatomical *.img file');
        anat_file= strcat(path,name)
    end
    
	% figure out name for anatomical file
    sz = size(anat_file);
    imgname = strcat(anat_file(1,1:sz(2)-4) , '.img');
    hdrname = strcat(anat_file(1,1:sz(2)-4) , '.hdr');
    %label=sprintf('%s\nANATOMY: ....%s', label, imgname(end-15:end));
	
    hdr = read_hdr(hdrname);
    anat_data = read_img2(hdr,imgname);
   
    % interpolate the paramter map to fit the anatomical one
    [x,y,z] = meshgrid(1:spm_hdr.ydim , 1:spm_hdr.xdim, 1:spm_hdr.zdim);
    [xi,yi, zi] = meshgrid(1:hdr.ydim , 1:hdr.xdim, 1:hdr.zdim);

    yi = yi * spm_hdr.ysize/hdr.ysize;
    xi = xi * spm_hdr.xsize/hdr.xsize;
    zi = zi * spm_hdr.zsize/hdr.zsize;

    spm_data2 = interp3(x,y,z, spm_data, xi,yi,zi,'nearest');
    spm_data2(find(isnan(spm_data2)))=0;
    %spm_data = spm_data2;

    %whos
    
    % threshold the map at the 50%
    if nargin ==0
        threshold = mean(mean(mean(spm_data2))) * 3
    end
    spm_data2(find(spm_data2 <= threshold )) = 0;
    
    % scale the maps to use the whole colormap
    anat_data = anat_data * 2^7 / max(max(max(anat_data)));
    spm_data2 = spm_data2 * 2^7 / max(max(max(spm_data2)));
    
    out_data = anat_data;
    out_data(find(spm_data2)) =  2^7 + spm_data2(find(spm_data2));
    
    %configure the colormap:
	mygray = [0:1/127:1]' * [1 1 1]; 
    
	myhot = [0:3:125]' * [1 0 0] ;
	tmp =   [0:3:125]' * [0 1 0] ;
	myhot = [myhot; tmp];
	tmp =   [0:3:125]' * [0 0 1];
	myhot =  [myhot;  tmp]/128;

	myhot(round(128/3): 128, 1) = 1;
	myhot(round(128*2/3):128,2) = 1;

 	mymap = [mygray; myhot];
    
	colormap(mymap)
	axis off
    
        
    d=out_data;
          
    if nargin < 5
        % determine which files contain the time series..:
        f='dummmyname';
        file=[];
        path = [];
        op=pwd;
        while f ~= 0
            [f p] = uigetfile('*.img','Select a file from the next run.  cancel when finished ');
            if f~=0
                file = [file ; f];
                path = [path ; p];
                cd(p);
            end
        end
        cd(op);
        numRuns = size(file,1);
    end
    
    
    xs = ceil(cx*spm_hdr.xdim/hdr.xdim);
    ys = ceil(cy*spm_hdr.ydim/hdr.ydim);
    zs = ceil(cz*spm_hdr.zdim/hdr.zdim);
    
    str=sprintf('(x,y,z)=  (%3.2f %3.2f %3.2f) mm,  \n (%3d %3d %3d)vx , \n val= %d',...
        hdr.xsize*(xs-spm_hdr.origin(1)), ...
        hdr.ysize*(ys-spm_hdr.origin(2)), ...
        hdr.zsize*(zs-spm_hdr.origin(3)), ...
        xs,ys,zs, ...
        spm_data(xs,ys,zs)*spm_scale );
    
    colordef black
    [ cx cy cz] 
    [fig1, fig2, fig3] =  ov(hdr,d,cx,cy,cz,0);
    subplot(221), title(str) %, xlabel(label)
    
    subplot(224), cla
    %if button==3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find a list of voxels within the sphere VOI

    roi_xyz = load(voxelFile);          

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make sure there are voxels above threshold
    if size(roi_xyz) ~= [0 0]  
        
        if nargin > 5
            numRuns=1;
            path=tseries_path;
            file=func_root;
        end
        % make a gaussian filter here	
        g = make_gaussian(50,4,100);
        
        tdata=[];
        for n=1:numRuns
            tmp = timeplot4(path(n,:), file(n,:), roi_xyz);
            
            % filtering stuff:
            % tmp = smoothdata2(tmp, TR, 0.009, 0.3, 3); 
            tmp = mydetrend(tmp);
            % load the mean (DC term) and divide by it
            % (compute percentage)
            load coeffs.mat
            mm = coeffs(end);
            tmp = 100 * tmp / mm;
            
            tmp = conv(g,tmp);
            tmp = tmp(50:end-50);
            tdata = [tdata tmp];
        end
        
        [ev_avg ev_std] = event_avg(tdata,onsets,window,10);
        
        subplot(221), hold on, plot(roi_xyz(:,2), roi_xyz(:,1),'g*');
        subplot(222), hold on, plot(roi_xyz(:,1), roi_xyz(:,3),'g*');
        subplot(223), hold on, plot(roi_xyz(:,2), roi_xyz(:,3),'g*');
        
        subplot 224, plot(ev_avg,'r');axis tight ; hold off; 
        subplot(224), title (sprintf('ROI at %d %d %d (%d vox)', xs , ys, zs, size(roi_xyz(2))))
        set(gca, 'Xtick',[0:2:window])
        
        
        tmp = [ev_avg ev_std];
        save avg.dat tmp -ASCII
        save tdata.dat tdata -ASCII
        save voxels.dat roi_xyz -ASCII
    end
    return
