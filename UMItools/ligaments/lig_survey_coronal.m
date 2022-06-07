function [h,M]=lig_survey_coronal(varargin)
% lig_survey_coronal_coronal(filename,varargin)

if ischar(varargin{1})
    P=load(varargin{1});
end

if isstruct(varargin{1})
    P=varargin{1};
end


L=P.Lagrange;

%%
hf=0;
mask=P.mask_hr.*P.mask; %and the masks of the GE3D data and the DESTE data

%identify a default center of the sample, using the mask;
si=size(mask);
mc1=round(si(1)/2);
testslice=squeeze(mask(mc1,:,:));
mc2=round(sum(testslice,2)'*(1:si(2))'/sum(testslice(:)));
mc3=round((1:si(3))*sum(testslice,1)'/sum(testslice(:)));
mask(mask==0)=NaN;
cmap=jet(256);

if any(strcmp(varargin,'ind'))
    dummy=find(strcmp(varargin,'ind'));
    slivec=varargin{dummy+1};
elseif nargin>1 && isnumeric(varargin{2});
    slivec=varargin{2};
else
    slivec=mc3+(-21:1:21);
end

if any(strcmp(varargin,'slim'))
    dummy=find(strcmp(varargin,'slim'));
    slim=varargin{dummy+1};
else
    slim=[-0.05 0.05];
end
    

if isfield(P.HIRES,'orig_data')
	[signal_estimate,noise_estimate]=estimate_noiselevel(P.HIRES.orig_data.image);
    clim=[noise_estimate 2.5*signal_estimate];
else
	[signal_estimate,noise_estimate]=estimate_noiselevel(P.HIRES.magnitude);
    clim=[noise_estimate 2.5*signal_estimate];
end

%% go through the slivecs
for sli=slivec
    
	% get the high resolution data axes
	if isfield(P.HIRES,'orig_data')				%if the HIRES data has better resolution
		ax1=P.HIRES.orig_data.pars.axis2;
		ax2=P.HIRES.orig_data.pars.axis1;
        picbw=(abs(P.HIRES.orig_data.image(:,:,sli,1)));
        
        if isfield(P.HIRES.orig_data,'mask' & any(strcmp(varargin,'mask')))
            picbw=picbw.*squeeze(P.HIRES.orig_data.mask(:,:,sli,1));
        end
		%clim=[0 20000];
	else
		ax1=P.HIRES.axis2;
		ax2=P.HIRES.axis1;
		picbw=abs(P.HIRES.magnitude(:,:,sli,1));
		%clim=[0 20000];
	end
		
	
    hf=hf+1;
    h(hf)=figure('position',[10 10 1650 950]);
    colormap(gray);
    set(gcf,'PaperPositionMode','manual','PaperOrientation','landscape');
    for enc=1:3
        
        tsubplot(3,7,(enc-1)*7 +1+enc,10);
        imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
        axis image
        %grid on
        set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
        title(['L ' num2str(enc) num2str(enc)],'color','k');
        
        %overlay the diagonal Lagrange strains
        lig_heatoverlay(L(:,:,sli,enc,enc),'clim',slim,'colorbar','noise',noise_estimate);
    end
    
    % L12
    tsubplot(3,7,3,10);
    imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
    axis image
    %grid on
    set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
    title(['L 12'],'color','k');
    %overlay the Lagrange strains
    lig_heatoverlay(L(:,:,sli,1,2),'clim',slim,'colorbar','noise',noise_estimate);
    
    % L13
    tsubplot(3,7,4,10);
    imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
    axis image
    %grid on
    set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
    title(['L 13'],'color','k');
    %overlay the Lagrange strains
    lig_heatoverlay(L(:,:,sli,1,3),'clim',slim,'colorbar','noise',noise_estimate);

    % L23
    tsubplot(3,7,11,10);
    imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
    axis image
    %grid on
    set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
    title(['L 23'],'color','k');
    %overlay the Lagrange strains
    lig_heatoverlay(L(:,:,sli,2,3),'clim',slim,'colorbar','noise',noise_estimate);
   
    %plot Q, calculated by derivs to strain
    tsubplot(2,5,4,4);
    imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
    axis image;
    title('Q equiv strain')
    %set(gca,'position',[0.1   0.0267    0.1680    0.2800]);
    set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
    lig_heatoverlay(P.Q(:,:,sli),'clim',slim*2,'colorbar','noise',noise_estimate);
    
    
    cutvec=(12:-1:1)*15;
    %plot V, calculated by derivs to strain
    tsubplot(2,5,9,4);
    imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
    axis image;
    title('Volm. strain');
    %set(gca,'position',[0.2    0.0267    0.1680    0.2800]);
    
    %set(gca, 'XColor', 'w','Ycolor','k','gridlinestyle','-','gridcolor',[1 1 1],'gridalpha',0.5,'YTick',round(10*ax2(cutvec))/10,'XTick',[],'Xticklabel','', ...
    %   'Tickdir','out');
    set(gca, 'XColor', 'w','Ycolor','k','gridlinestyle','-','YTick',round(10*ax2(cutvec))/10,'XTick',[],'Xticklabel','', ...
       'Tickdir','out');
    
    lig_heatoverlay(P.V(:,:,sli),'clim',slim*2,'colorbar','noise',noise_estimate);
    axis on;
    grid on;
    
    % angle plots in the left column
	titlevec={'ud'; 'lr'; 'io'};
	for enc=1:3
		tsubplot(3,6,(enc-1)*6+1,10);
		imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
		axis image
		%grid on
		set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
		title([titlevec{enc} ',\lambda=' num2str(P.lambda(enc))]);
		
		%overlay the diagonal Lagrange strains
		lig_angleoverlay(P.data(:,:,sli,enc),'colorbar','noise',noise_estimate);
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Axial images with a line indicating the slice
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% get the high resolution data axes
	if isfield(P.HIRES,'orig_data')				%if the HIRES data has better resolution structure (sems)
		ax1=P.HIRES.orig_data.pars.axis2;
		ax2=P.HIRES.orig_data.pars.axis3;
		ax3=P.HIRES.orig_data.pars.axis1;
		picbw=abs(P.HIRES.orig_data.image);
		picbw=permute(picbw,[3 2 1]);
		%clim=[0 20000];
	else
		ax1=P.HIRES.axis2;
		ax2=P.HIRES.axis3;
		ax3=P.HIRES.axis1;
		picbw=abs(P.HIRES.magnitude);
		picbw=permute(picbw,[3 2 1]);
		%clim=[0 20000];
    end
    
    %plot axial cuts, with a line indicating the slice
    ny=ceil(numel(cutvec))/2;
    for jj=1:numel(cutvec)/2
        tsubplot(ny,10,jj*10-1,15);
        mc1=cutvec(2*(jj-1)+1);
        % plot the magnitudes in the 5th column
        %MM=permute(squeeze((abs(P.HIRES.orig_data.image(mc1,:,:)))),[2 1]);
		MM=picbw(:,:,mc1);
        imagesc(ax1,ax2,MM,[clim(1) clim(2)]);
        axis image
        set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','YDir','reverse','Yticklabel','','Xticklabel','');
        hold on;
        plot([ax1(1) ax1(end)],ax2(sli)*[1 1],'w');
        title(num2str(round(10*ax3(mc1))/10));
        
        tsubplot(ny,10,jj*10,15);
        mc1=cutvec(2*jj);
        % plot the magnitudes in the 5th column
        MM=picbw(:,:,mc1);
        imagesc(ax1,ax2,MM,[clim(1) clim(2)]);
        axis image
        set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','YDir','reverse','Yticklabel','','Xticklabel','');
        hold on;
        plot([ax1(1) ax1(end)],ax2(sli)*[1 1],'w');
        title(num2str(round(10*ax3(mc1))/10));
	end
	
	
	
    
    set(gcf,'Renderer','OpenGL');
    
    M(hf)=getframe(gcf);
    
end


%save('ligmovie.mat','M');
%for k=64:-1:1, mx(:,:,k,1) = M(k).cdata(:,:,1); mx(:,:,k,2) = M(k).cdata(:,:,2); mx(:,:,k,3) = M(k).cdata(:,:,3); end


if any(strcmp(varargin,'pdf'))
		filename='DESTE';
		for enc=1:3
			filename=[filename '_' P.timest{enc}];
		end
		
		filename=[filename '_wiggle=' num2str(mean(P.wiggle)) 'mm'];
		
		filename=[filename '.pdf'];
		
		
        pdfappend(h,filename,'size',[10 8]);
    end
end
