function [im, h]= lightbox(root, wscale, rows, imgNum)
%function im = lightbox(root, [wscale],[ rows],[imageNum])
%
%   (c) 2005 Luis Hernandez-Garcia
%   University of Michigan
%   report bugs to:  hernan@umich.edu
%
% This program displays the slices in a 3D data set OR a time series
% root :  either the name of a file in a time series (could be asingle file)
% 	or could also be a 3D matrix that you want to display in slices
%	the program checks to see if it's a string or a matrix...
% wscale: window scale factor for display
% rows: number of rows of slices in the lightbox
%

global RT_MODE

if nargin<4
    num=1;
else
    num=imgNum;
end

if isstr(root)
    [im,h] = read_img(root);
    if isfield(h,'magic')
        h=nii2avw_hdr(h);
    end
    
    if h.tdim == 1
        im = reshape(im,h.xdim, h.ydim, h.zdim);
    else
        im = reshape(im(num,:) ,h.xdim, h.ydim, h.zdim);
    end
else
    im = root;
    h.xdim=size(im,1);
    h.ydim=size(im,2);
    h.zdim=size(im,3);
    
end

if nargin==1
    rows=[];
    wscale=[];
end
if isempty(rows)
    rows = floor(sqrt(h.zdim)) ;
end
M=[];
cols=ceil(h.zdim/rows);

for r=1:rows
    Mrow = [];
    
    for c=1:cols
        sl = c + cols*(r-1);
        if sl<=h.zdim
            Mrow = [Mrow  im(:,:,sl)'];
        else
            Mrow = [Mrow  zeros(h.ydim, h.xdim)];
        end
        
    end
    M = [M ; Mrow];
end

imagesc(M);

axis image
axis xy

if isempty(RT_MODE) || (RT_MODE ~=1)
    colormap(gray(256));
    grid off
    axis off
else
    set(gca,'Position',[0.01 0.01 0.99 0.99])
end

if ~isempty(wscale)
    caxis(wscale);
end

colorbar

if isstr(root)
    set(gcf,'Name',root)
end
return

function out = show(a,windfact)
% usage .. show(a,f);
% displays matrix "a" as a greyscale image in the matlab window
% and "f" is the optional window factors [min,max]

global RT_MODE

if isempty(RT_MODE) || (RT_MODE ~=1)
    doColorBar=1;
else
    doColorBar=0;
end

if exist('windfact') == 0,
    amin = min(min(a));
    amax = max(max(a));
    minmax = [amin,amax];
    a = (a  - amin);
else
    amin = windfact(1);
    amax = windfact(2);
    minmax = [amin,amax];
    a = (a  - amin);
    
end


%a(a>amax) = amax;


a = a./(amax-amin).*255;
a(a < 0)=0;
a(a>255) = 255;

imageHandle = image(a);


if doColorBar
    ncbar_labels=5;
    c1 = colorbar;
    set(c1,'YTick',linspace(1,256,ncbar_labels),...
        'YTickLabel',linspace(amin,amax,ncbar_labels),...
        'FontSize',12);
end

if nargout==0 && isempty(RT_MODE),
    disp(['min/max= ',num2str(minmax(1)),' / ',num2str(minmax(2))]);
end;
