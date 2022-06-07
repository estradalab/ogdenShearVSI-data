function out=lig_compare_ge_se(varargin)
% function to read in processed se and ge echo dat for stretching
% experiments, plotting geometries side by side
% inputs are time stamps

for jj=1:nargin;
    %make the time stamp string
    if isnumeric(varargin{jj})
        timechar{jj}=['0000' num2str(varargin{jj})];
        timechar{jj}=timechar{jj}(end-3:end);
    else
        timechar{jj}=varargin{jj};
    end
        
    
    D{jj}=load(['HIRES_' timechar{jj} '_magnitude.mat']);
end

figure('position',[600 600 400*nargin 400]);

for kk=1:2:192;
    
    for jj=1:nargin;
        subplot(1,nargin,jj);
        
        pulsesequence=D{jj}.params.seqfil;
        
        switch pulsesequence
            case 'ge3d'
                axialslice=abs(squeeze(D{jj}.compleximage(kk,:,:))');
                xaxis=D{jj}.axis2;
                yaxis=D{jj}.axis3;
                titlestr=[timechar{jj} ' ge3d'];
            case 'sems'
                axialslice=abs(squeeze(D{jj}.orig_data.image(kk,:,:)))';
                xaxis=D{jj}.orig_data.pars.axis2;
                yaxis=D{jj}.orig_data.pars.axis3;
                titlestr=[timechar{jj} ' sems'];
        end
        
        imagesc(xaxis,yaxis,axialslice);
        axis image;
        grid on;
        set(gca,'Ydir','normal','gridcolor',[ 1 1 1],'gridalpha',0.5, ...
            'xtick',[-7.5:2.5:7.5],'ytick',[-7.5:2.5:7.5]);

        title(titlestr);
    end
    pause(0.5)
end
        
        