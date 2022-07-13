function out=lig_heatoverlay(varargin)
% out=lig_heatoverlay(heatmatrix [,clim]);
% options:
% 'clim' defines upper and lower bound of heat, say [-pi pi]
% 'noise',noiselevel    : passes noiselevel of underlying (BW) image; --> opacity =0 where SNR~1
% 'flat'    uses the full color range, including where heat is zero 
%

set(gca,'colormap',gray(128));
Heat=varargin{1};
Heat(isnan(Heat(:)))=0;

% get colorlimits from range of heat data, or varargin
if any(strcmp(varargin,'clim'))
    ind=find(strcmp(varargin,'clim'));
    clim=varargin{ind+1};
else
    clim=[min(Heat(:)) max(Heat(:))];
end


%get the image data contained in the current axis
children=get(gca,'children');
% child #1 is the current BW image
ih=children(end); %the last child is the first image drawn to the axis

% read the image magnitude, get its signal and noise estimates
I1.CData=abs(get(ih,'CData'));  %these are the actual data values
I1.XData=get(ih,'XData');
if numel(I1.XData)==2
    I1.XData=I1.XData(1):I1.XData(2);
end
I1.YData=get(ih,'YData');
if numel(I1.XData)==2
    I1.YData=I1.YData(1):I1.YData(2);
end

% get noise from BW image or varargin input
if any(strcmp('noise',varargin))
    ind=find(strcmp('noise',varargin));
    noiseA=varargin{ind+1};
else
    [~,noiseA]=estimate_noiselevel(I1.CData);
end
    

%get the size of the image to be overlaid, and make a coordinate grid
siA=size(I1.CData);
axA1=(1:siA(1))/siA(1); axA2=(1:siA(2))/siA(2);
[AG1,AG2]=ndgrid(axA1,axA2);

% get the size and make temp axis for the overlay handed down
Heat=real(Heat);
siHeat=size(Heat);
axHeat1=(1:siHeat(1))/siHeat(1); axHeat2=(1:siHeat(2))/siHeat(2);
[HeatG1,HeatG2]=ndgrid(axHeat1,axHeat2);

%regrid the overlay
iHeat=interpn(HeatG1,HeatG2,Heat,AG1,AG2);
iHeat(isnan(iHeat(:)))=0;

%make a color image based on the heat map
col_overlay=generateRGB(iHeat,clim);

%define a weighting function for the color  overlay
% it may depend on 
% (1) the noise level of the BW image
% (2) the signal level of BW image, or the log signal level
% (3) the values in the heatmap 

% first, define amplitude weighting as a function of SNR in the BW picture
% 0 for SNR<1, 1 for SNR>3, logarithmic in between; rejects empty space
amplitudeweighting=min(1,max(0.05,log10(I1.CData/noiseA)));

% second, define heat (STRAIN) weighting, to reject coloration where there is no heat
% 0 for Heat less than 10% of max range, 1 for Heat > 27% of max range 
heatrange=max(abs(clim));
if ~any(strcmp('flat', varargin))
    heatweighting=ones(size(amplitudeweighting));
else
    heatweighting=min(1,max(0.2,log10(abs(iHeat/(heatrange/40)))));
end


% range fro 0 to 1
weighting=heatweighting.*amplitudeweighting;  %the color overlays get weighted by image intensity and by the sign and magnitude of the heat map

col_opacity=abs(weighting/2);  %factor 2 guarantees at most 50% opacity


%get the current position
pos=get(gca,'Position');
hold on;
%overlay the color map, give it opacity proportional to the weighting
%factor
hc=imshow(col_overlay);
set(hc,'XData',I1.XData);
set(hc,'YData',I1.YData);
set(hc,'alphadata', col_opacity);
hold off;

if any(strcmp(varargin,'colorbar'))
    currentaxes=gca;
    hcb=colorbar;  %axis handle of the colorbar
    pos=get(hcb,'Position');
    set(hcb,'visible','off');
    axes('Position',pos,'XTicklabel','');
    dummy=(1:128)'*[1 1 1 1 1];
    cm=generateRGB(dummy);
    
    cb=imshow(cm);
    axis on;
    ytick=0.5:32:128.5;
    yticklabel=min(clim)+[0:4]*abs(clim(2)-clim(1))/4;
    
     set(gca,'YDir','normal','Yaxislocation','right','XTicklabel','','XTick',[],...
        'YTick',ytick,'yticklabel',yticklabel,'Tickdir','in');
    axes(currentaxes);
end

if any(strcmp(varargin,'diagnostic'))
	gca,
    
	figure;
	h=subplot(2,5,1);
	imagesc(ih.XData,ih.YData,ih.CData);
	axis image;
	colormap gray;
	title('anatomy');
	
	subplot(2,5,2);
	dummy=sqrt(max(0,log10(I1.CData/noiseA)));
	imagesc(ih.XData,ih.YData,dummy);
	axis image;
	ht=title('wt=[ln(SNR))]^{0.5}');
	
	subplot(2,5,7);
	imagesc(ih.XData,ih.YData,col,[-0.3 0.3]);
	axis image;
	title('equiv. strain Q');
	
	subplot(1,4,3);
	imagesc(ih.XData,ih.YData,weighting);
	axis image;
	title('Q \times wt');

	subplot(1,4,4);
	imagesc(ih.XData,ih.YData,I1.CData);
	axis image;
	hold on;
	hc=imshow(col);
	set(hc,'XData',I1.XData);
	set(hc,'YData',I1.YData);
	set(hc,'alphadata', col_opacity);
	hold off;
	title('overlay');
end

out=gcf;

function out=generateRGB(varargin)
%takes the inmatrix MxN and turns it into a MxNx3 RGB map
%for imshow overlay on the current map
inmatrix=varargin{1};
if nargin>1
    %clim defined
    clim=varargin{2};
    mincl=min(varargin{2});
    maxcl=max(varargin{2});
    inmatrix(inmatrix<mincl)=mincl;
    inmatrix(inmatrix>maxcl)=maxcl;
else
    clim=[min(inmatrix(:)) max(inmatrix(:))];
end

load spectral;

if clim(1)==0; %ascending colormap from half point
    cmap=spectral(128:255,:);
elseif clim(2)==0; %ascending colormap from bottom to half point;
    cmap=spectral(1:128,:);
else
    cmap=spectral;
end

C = cmap;  % Get the figure's colormap.

% 
% cm=jet(128);
% 
% C = jet(128);  % Get the figure's colormap.
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(clim(1),clim(2),L),1:L,inmatrix));
out= reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.

