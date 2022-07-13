function result = orthomov(root )
% function result = orthomov(root)
%
% this function thakes a file in analyze format and 
% displays orthogonal sections thru the planes interscting
% at x,y,z
global wscale

origFig = figure;
movieFig = figure;
figure(origFig);
colormap(gray);
roi=0;

% configure colormap
my_map=(0:255)';
my_map=[my_map my_map my_map]/256;
colormap(my_map);

% determine which files to use
if nargin==0
    [fname path] = uigetfile('*.img','Select time series 4D file');
    root= strcat(path,fname);
end

% display the first one and let user determine planes
[whole h] = read_img(root);
numImgs = h.tdim;


x=round(h.xdim/2);
y=round(h.ydim/2);
z=round(h.zdim/2);

d = whole(1,:);
% scale image to fit colormap
lo = min(min(min(d)));
hi = max(max(max(d)));
if ~isempty(wscale)
    hi=wscale(2);
    lo = wscale(1);
end
range = hi-lo;
dd = (d-lo)*256/range;

i=1;
% interactive display ...
while i >= -10
   figure(movieFig)
   colormap(my_map)

   d = whole(2,:);
   d = (d-lo)*256/range;
   d = reshape(d,h.xdim, h.ydim, h.zdim);
   [fig1, fig2, fig3] = ov(h,d,x,y,z, roi);
    
    [i j button] = ginput(1)
    
    i=round(i);
    j=round(j);
    fig = floor(gca);
    
    if button == 1
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
    end
        
        
        %play movie of selected view
        if button == 3

                title('saving ortho.avi')
        	aviF = avifile('ortho.avi')
                aviF.fps = 4;
        	clear M;
        end
        for idx = 2:numImgs
            d = whole(idx,:);
            d = (d-lo)*256/range;
            d = reshape(d,h.xdim, h.ydim, h.zdim);
            [fig1, fig2, fig3] = ov(h,d,x,y,z, roi);
	    drawnow
	    if button==3
                  M(idx) = getframe(gcf);
                  aviF=addframe(aviF, M(idx));
                  fprintf('Frame : %d  ',idx);
            end
        end
	if button==3
             aviF=close(aviF);
        end
        %fps = input('Enter Frame Rate (in FPS): ');
        fps = 4	;
        figure(origFig);
end

return
