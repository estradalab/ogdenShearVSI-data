function [outfit,roicell,synthimage]=brain_mefit(varargin)
%spine multi echo fit
%[outfit,roicell,synthimage]=brain_mefit(structure_returned_by_varianms [, sliceindex-for_picking_ROI] [,roi] [,'oddeven']);

data=varargin{1};
imagedata=abs(data.image);

xaxis=data.pars.xaxis;
yaxis=data.pars.yaxis;
sli=varargin{2};                %slice index

%loop to check what has been passed down
for jj=1:nargin;
    structflag(jj)  =isstruct(varargin{jj});
    cellflag(jj)    =iscell(varargin{jj});
    realflag(jj)    =isreal((varargin{jj}));
    %optional other type determinats
end


if any(cellflag);
    dummy=find(cellflag);
    for jj=1:numel(dummy);
        testcell=varargin{dummy(jj)};
        if strcmp(testcell{1},'ROI limits');
            xlim=testcell{2};
            ylim=testcell{3};
            xwinvec=testcell{4};
            ywinvec=testcell{5};
            
            roicell=testcell;
        end
    end
else
    [digiout, xlim, ylim]=pickroi(imagedata(:,:,sli,1));
    xwinvec=xlim(1):xlim(2);
    ywinvec=ylim(1):ylim(2);
    roicell={'ROI limits'; xlim; ylim;xwinvec;ywinvec};
end

%Parameter setup mems or mgems
datatofit=squeeze(imagedata(ywinvec,xwinvec,:,:));  %default is fit all the data
if (strcmp(data.pars.seqfil,'mgems') || strcmp(data.pars.seqfil,'mge3d')) ;
    si=size(datatofit);
    taxis=data.pars.te+[0:data.pars.ne-1]*data.pars.te2;
    
    %restrict the times over which to fit
    restrictT2starflag=false;
    if restrictT2starflag;
        tindvec=(taxis<25e-3);
        taxis=taxis(tindvec);
        datatofit=squeeze(imagedata(ywinvec,xwinvec,:,tindvec));
    end
   
    
    figuretitle=['T2* decay'];
    
    TREF=0.03;
    taulim=[0 0.8];
    TEvec=[5 10 15]/1000;
    
elseif strfind(data.pars.seqfil,'mems');
    si=size(datatofit);
    taxis=data.pars.te*[1:data.pars.ne];
    
    %restrict the times over which to fit
    restrictT2flag=true;
    if restrictT2flag;
        tindvec=(taxis<300e-3);
        taxis=taxis(tindvec);
        datatofit=squeeze(imagedata(ywinvec,xwinvec,:,tindvec));
    end

    
    figuretitle=['MEMS T2 decay'];
    taulim=[0 0.2];
    TREF=0.1; 
    TEvec=[50 100 200]/1000;
    
end

[nofRO, nofPE,nofSL,nofTimes]=size(datatofit);
amplitude=zeros(nofRO,nofPE,nofSL);
tau=zeros(nofRO,nofPE,nofSL);
synthimage=zeros(nofRO,nofPE,nofSL,numel(TEvec));
%fit the data to exponential decay
if any(strcmp(varargin,'oddeven'));
    for jj=1:nofSL;
        outfit=matrix_expfit(squeeze(datatofit(:,:,jj,:)),taxis, 'oddeven');
        amplitude(:,:,jj)=outfit.amplitude;
        tau(:,:,jj)=outfit.Tau;
        
        %calculate the synthetic image
        for te=1:numel(TEvec);
            synthimage(:,:,jj,te)=outfit.amplitude.*exp(-TEvec(te)./outfit.Tau);
        end 
    end
else
    for jj=1:nofSL;
        outfit=matrix_expfit(squeeze(datatofit(:,:,jj,:)),taxis);
        amplitude(:,:,jj)=outfit.amplitude;
        tau(:,:,jj)=outfit.Tau;
        
        % calculate the synthetic image
        for te=1:numel(TEvec);
            synthimage(:,:,jj,te)=outfit.amplitude.*exp(-TEvec(te)./outfit.Tau);
        end
    end
end

%characterize the distribution of pixel intensities in all images
for sl=1:nofSL;
    test=amplitude(:,:,sl);
    meanintensity(sl)=mean(test(:));
    stdintensity(sl)=std(test(:));
end

cutoff=mean(meanintensity-stdintensity)/2;
tau(amplitude<cutoff)=0;
tau(tau<0)=0;
tau(tau>5)=0;


%if plotflag
%plot amplitdue, Tau, and calculated synthetic images, for all slice
% plot figures 

plotflag=false;
if any(strcmp(varargin,'plot'))
    plotflag=true;
end

if plotflag;
    %amplitude plot
    clim=[0 2.5*mean(amplitude(:))];
    plotms(amplitude,'clim',clim);
    colormap gray(256);
    set(gcf,'Name',[figuretitle ', amplitude']);
    
    %tau plot
    plotms(tau,'clim',taulim');
    colormap gray(256);
    set(gcf,'Name',[figuretitle ', tau']);
end

clear outfit
outfit.amplitude=amplitude;
outfit.Tau=tau;
