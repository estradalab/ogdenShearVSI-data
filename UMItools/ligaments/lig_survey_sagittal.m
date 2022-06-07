function h=lig_survey_sagittal(varargin)
% lig_survey_coronal_sagittal(filename,varargin)

if ischar(varargin{1})
    P=load(varargin{1});
end

if isstruct(varargin{1})
    P=varargin{1};
end


L=P.Lagrange;

% %%
hf=0;
if isfield(P,'mask_hr');
    mask=P.mask_hr.*P.mask; %and the masks of the GE3D data and the DESTE data
else
    mask=P.mask;
end
%mask=P.mask;

%Default center of the sample, using the mask;
si=size(mask);
mm=round(si/2);  mc1=mm(1); mc2=mm(2); mc3=mm(3);

if any(strcmp(varargin,'ind'))
    dummy=find(strcmp(varargin,'ind'));
    slivec=varargin{dummy+1};
else
    slivec=mc2+(-21:3:21);
end

if isfield(P.HIRES,'orig_data')
	[signal_estimate,noise_estimate]=estimate_noiselevel(P.HIRES.orig_data.image);
    clim=[noise_estimate 5*signal_estimate];
else
	[signal_estimate,noise_estimate]=estimate_noiselevel(P.HIRES.magnitude);
    clim=[noise_estimate 5*signal_estimate];
end

%% go through the slivecs
for sli=slivec
    
	% get the high resolution data axes
	if isfield(P.HIRES,'orig_data')				%if the HIRES data has better resolution
		ax1=P.HIRES.orig_data.pars.axis3;
		ax2=P.HIRES.orig_data.pars.axis1;
        %the slice index into the strain data is predicated on a data set
        %size of [192 64 64], reached by upsampling. The
        %HIRES.orig_data.image may have better resolution int the phase
        %encode direction (readout and slice are forced to 192 and 64 in lig_deste_3d) 
        %so the slice index for sagittal slices must be adapted accordingly
        si1=size(P.data);
        si2=size(P.HIRES.orig_data.image);
		picbw=squeeze(abs(P.HIRES.orig_data.image(:,round(sli*si2(2)/si1(2)),:,1)));
        if isfield(P.HIRES.orig_data,'mask')
             picbw=picbw.*squeeze(P.HIRES.orig_data.mask(:,round(sli*si2(2)/si1(2)),:,1));
        end
		%clim=[0 20000];
    else
        
		ax1=P.HIRES.axis3;
		ax2=P.HIRES.axis1;
		picbw=squeeze(mask(:,sli,:).*(abs(P.HIRES.magnitude(:,sli,:,1))));
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
        lig_heatoverlay(squeeze(L(:,sli,:,enc,enc)),'clim',[-0.1 0.1],'colorbar','noise',noise_estimate);
    end
    
    % L12
    tsubplot(3,7,3,10);
    imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
    axis image
    %grid on
    set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
    title(['L 12'],'color','k');
    %overlay the Lagrange strains
    lig_heatoverlay(squeeze(L(:,sli,:,1,2)),'clim',[-0.1 0.1],'colorbar','noise',noise_estimate);
    
    % L13
    tsubplot(3,7,4,10);
    imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
    axis image
    %grid on
    set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
    title(['L 13'],'color','k');
    %overlay the Lagrange strains
    lig_heatoverlay(squeeze(L(:,sli,:,1,3)),'clim',[-0.1 0.1],'colorbar','noise',noise_estimate);

    % L23
    tsubplot(3,7,11,10);
    imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
    axis image
    %grid on
    set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
    title(['L 23'],'color','k');
    %overlay the Lagrange strains
    lig_heatoverlay(squeeze(L(:,sli,:,2,3)),'clim',[-0.1 0.1],'colorbar','noise',noise_estimate);
   
    %plot Q, calculated by derivs to strain
    tsubplot(2,5,4,4);
    imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
    axis image;
    title('Q equiv strain')
    %set(gca,'position',[0.1   0.0267    0.1680    0.2800]);
    set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
    lig_heatoverlay(squeeze(P.Q(:,sli,:)),'clim',[-0.2 0.2],'colorbar','noise',noise_estimate);
    
    
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
    
    lig_heatoverlay(squeeze(P.V(:,sli,:)),'clim',[-0.2 0.2],'colorbar','noise',noise_estimate);
    axis on;
    grid on;
    
    %plot the three phase data plots in color, using the HIRES.magnitude as
	%the background, i.e. the BW image with the same resolution as the
	%angle image
    titlevec={'ud'; 'io'; 'lr'};
	for enc=1:3
		tsubplot(3,6,(enc-1)*6+1,10);
		imagesc(ax1,ax2,picbw,[clim(1) clim(2)]);
		axis image
		%grid on
		set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
		title([titlevec{enc} ',\lambda=' num2str(P.lambda(enc))]);
		
		%overlay the phases
		lig_angleoverlay(squeeze(angle(P.data(:,sli,:,enc))),'colorbar','noise',noise_estimate);
	end
    
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Axial images with a line indicating the slice
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% get the high resolution data axes
	if isfield(P.HIRES,'orig_data')				%if the HIRES data has better resolution structure (sems)
		ax1=P.HIRES.orig_data.pars.axis3;
		ax2=P.HIRES.orig_data.pars.axis2;
		ax3=P.HIRES.orig_data.pars.axis1;
		picbw=abs(P.HIRES.orig_data.image);
		picbw=permute(picbw,[3 2 1]);           %%%%%
		%clim=[0 20000];
	else
		ax1=P.HIRES.axis3;
		ax2=P.HIRES.axis2;
		ax3=P.HIRES.axis1;
		picbw=abs(P.HIRES.magnitude);
		picbw=permute(picbw,[3 2 1]);           %%%%%
		%clim=[0 20000];
	end
	
	
	
	
    
    %plot axial cuts, with a line indicating the slice
    ny=ceil(numel(cutvec))/2;
    for jj=1:numel(cutvec)/2
        
        %left panel
        tsubplot(ny,10,jj*10-1,15);
        mc1=cutvec(2*(jj-1)+1);
        % plot the magnitudes in the 5th column, 
        %picbw is permuted to [3 2 1] just above
		MM=picbw(:,:,mc1);
        imagesc(ax1,ax2,MM,[clim(1) clim(2)]);
        axis image
        set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','YDir','reverse','Yticklabel','','Xticklabel','');
        hold on;
        plot(ax1(sli)*[1 1],[ax2(1) ax2(end)],'w');
        title(num2str(round(10*ax3(mc1))/10));
        
        %right panel
        tsubplot(ny,10,jj*10,15);
        mc1=cutvec(2*jj);
        % plot the magnitudes in the 5th column
        MM=picbw(:,:,mc1);
        imagesc(ax1,ax2,MM,[clim(1) clim(2)]);
        axis image
        set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','YDir','reverse','Yticklabel','','Xticklabel','');
        hold on;
        plot(ax1(sli)*[1 1],[ax2(1) ax2(end)],'w');
        title(num2str(round(10*ax3(mc1))/10));
	end
	
    
%     M(sli)=getframe(gcf);
    
end