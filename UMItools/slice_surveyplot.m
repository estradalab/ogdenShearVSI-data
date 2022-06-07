function h=slice_surveyplot(varargin)
pic=varargin{1}.image;

if any(strcmp(varargin,'col'))
    cmap=jet(64);
else
    cmap=gray(64);
end
%auto choose the layout
si=size(pic);

if ~any(strcmp(varargin,'clim'))
    maxamp=abs(max(pic(:)));
    clim=[0 maxamp];
else
    ind=find(strcmp(varargin,'clim'));
    clim=varargin{ind+1};
end

switch ndims(pic)
    case 2
        nfig=1;
        nx=1;
        ny=1;
        nimperfigure=1;
    case 3						%stack of slices or time stamped single slice
        nfig=1;
        nx=ceil(sqrt(si(3)));
        ny=ceil(si(3)/nx);
        nimperfigure=si(3);	%number of images per figure
    case 4
        nfig=1;
        nx=si(3);
        ny=si(4);
        nimperfigure=nx*ny;
    case 5
        nfig=si(5);
        nx=si(3);
        ny=si(4);
        nimperfigure=nx*ny;
end

if nargin>1 && any(strcmp(varargin,'square'))
    sipic=size(pic);
    npanels=prod(sipic(3:end));
    nx=ceil(sqrt(npanels));
    ny=ceil(npanels/nx);
    nfig=1;
    nimperfigure=npanels;
end

if any(strcmp(varargin,'layout'))
    ind=find(strcmp(varargin,'layout'));
    layoutvec=varargin{ind+1};
    if numel(layoutvec)==2
        nfig=1;
        ny=layoutvec(1);
        nx=layoutvec(2);
        nimperfigure=nx*ny;
    elseif numel(layoutvec)==3
        ny=layoutvec(1);
        nx=layoutvec(2);
        nfig=layoutvec(3);
        nimperfigure=nx*ny;
    end
    
   %define an aspect ratio for the figure; as square won't do
   %the default plot is image without units, so square pixels
   aspectratio=(si(1)*ny)/(si(2)*nx);
   figsize=[600/aspectratio 600]; 
else
    figsize=[800 800];
end

if nargin>1 && any(strcmp(varargin,'slice')) %Slicewise
    slfl=find(strcmp(varargin,'slice'));
    if nargin>slfl && isnumeric(varargin{slfl+1}) %maxamplitude specified;
        maxamp=varargin{slfl+1};
    end
    sipic=size(pic);
    npanels=prod(sipic(3:end))/varargin{1}.pars.ns;
    nx=ceil(sqrt(npanels));
    ny=ceil(npanels/nx);
    nfig=varargin{1}.pars.ns;
    nimperfigure=npanels;
    for ff=1:nfig;
        h(ff)=figure('Name', ['slice ' num2str(ff)],'NumberTitle', 'off','Position',[ 150          100        figsize]);
        %h(ff)=figure('Position',[ 350          100        1000         1000]);
        colormap(cmap);
        for ii=1:nimperfigure;
            subplot(ny,nx,ii);
            imagesc(abs(pic(:,:,ff,ii)));
            axis image
            set(gca,'clim',clim,'XTick', [],'YTick',[]);
            title(['Echo #' num2str(ii)]);
        end
    end
else
    %plot it
    for ff=1:nfig
        %h(ff)=figure('Name', pc.comment,'NumberTitle', 'off');
        h(ff)=figure('Position',[ 150          100        figsize]);
        colormap(cmap);
        for ii=1:nimperfigure
            subplot(ny,nx,ii);
            imagesc(abs(pic(:,:,ii+(ff-1)*nimperfigure)));
            axis image
            set(gca,'clim',clim,'XTick', [],'YTick',[]);
            %title(['Image #' num2str(ii+(ff-1)*nimperfigure)]);
        end
    end
end
