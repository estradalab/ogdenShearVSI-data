function [meanVal, stdVal, indices] = getROIvals(roi, imfile, windfact,threshold)
% function [meanVal, stdVal, indices] = getROIvals(roi, imfile, windfact,threshold)

hdr = read_hdr(imfile);
D = read_img2(imfile);
glblMean = mean(D(find(D(:))));

x = round(hdr.xdim/2);
y = round(hdr.ydim/2);
z = round(hdr.zdim/2);

my_map=(0:255)';
my_map=[my_map my_map my_map]/256;
colormap(my_map);

D2 = D;
% scale image to fit colormap
if ~isempty('windfact'),
    dmin = min(D(:));
    dmax = max(D(:));
    D = (D  - dmin);
else
    dmin = windfact(1);
    dmax = windfact(2);
    D = (D  - dmin);
    D = D .* (D > 0);
end
D=D./(dmax-dmin).*256;


[fig1, fig2, fig3] =  ov(hdr,D,x,y,z,0);
[i j button] = ginput(1);

if strcmp(computer,'GLNX86')
    UPKEY = 30;
    DNKEY = 31;
    RTKEY = 29;
    LTKEY = 28;
end
if strcmp(computer,'MAC')
    UPKEY = 56;
    DNKEY = 50;
    RTKEY = 54;
    LTKEY = 52;
end

while button~=3  %3 is the right mouse button
    
    if(button==UPKEY) %30 = up arrow key to grow threshold
        roi=roi + 1
       
    elseif(button==DNKEY)  %31 = down arrow key to shrink threshold
        roi=roi-1
       
    elseif(button==RTKEY)  %29 = right arrow key to increase threshold by 10%
        threshold=threshold*1.1
       
    elseif(button==LTKEY)  %28 = left arrow key to decrease threshold by 10%
        threshold=threshold*0.9
       

    elseif(button==1)
        i=round(i);j=round(j);
        fig = floor(gca);
        switch(fig)
            case floor(fig1)
                x=j;
                y=i;
            case floor(fig2)
                z=j;
                x=i;
            case floor(fig3)
                y=i;
                z=j;
        end

        if exist('threshold') == 0,
            threshold=D2(x,y,z).*0.75
        end
    end
    xmin = max([1 x-roi]);
    xmax = min([hdr.xdim x+roi]);

    ymin = max([1 y-roi]);
    ymax = min([hdr.ydim y+roi]);

    zmin = max([1 z-roi]);
    zmax = min([hdr.zdim z+roi]);

    %nlist = (2*args.ROIsize + 1)^3; % total # of voxels within a cube
    [xx,yy,zz] = ndgrid([xmin:xmax], [ymin:ymax], [zmin:zmax]);
    nlist = length(xx(:));
    fx = reshape(xx,nlist,1);
    fy = reshape(yy,nlist,1);
    fz = reshape(zz,nlist,1);

    cube_xyz = [fx fy fz];   % a list of all voxels within the cube
    dist_cube = sqrt(...
        (hdr.xsize * (fx - x)).^2 + ...
        (hdr.ysize * (fy - y)).^2 + ...
        (hdr.zsize * (fz - z)).^2);    % a list of distances in mm

    sphere_xyz = cube_xyz(dist_cube <= roi, :); % a list of voxels within sphere
    num_voxels = size(sphere_xyz,1); % len_list is the # of voxels within sphere

    % keep only those above threshold in the spm_data
    func_xyz =  [];
    for v=1: num_voxels
        if D2(sphere_xyz(v,1), sphere_xyz(v,2), sphere_xyz(v,3)) ...
                >= threshold;
            func_xyz = [func_xyz; sphere_xyz(v,:)];
        end
    end

    if(sum(size(func_xyz))==0)
        meanVal = 0
        stdVal = 0

        [fig1, fig2, fig3] = ov(hdr,D,x,y,z,0);
    else
        vx = func_xyz(:,1);
        vy = func_xyz(:,2);
        vz = func_xyz(:,3);
        indices = sub2ind(size(D), vx, vy, vz);
        D3 = D2(:);

        meanVal = mean(D3(indices));
        stdVal = std(D3(indices));
         fprintf('\n radius = %0.2f mm , thresh = %0.2f, mean= %0.2f, STD = %0.2f ', ...
             roi, threshold, meanVal, stdVal);
        tmp = [meanVal stdVal glblMean];
        save voxels.dat func_xyz -ASCII
        save ROIstats.dat tmp -ASCII
        
        hold off
        [fig1, fig2, fig3] = ov(hdr,D,x,y,z,0);
        axes(fig3), hold on, plot(vy, vz,'gx');
        axes(fig2), hold on, plot(vx, vz,'gx');
        axes(fig1), hold on, plot(vy, vx,'gx');
        
        subplot(224), hist(D3(indices))
        title('ROI stats')
    end
    
    [i j button] = ginput(1);
end
