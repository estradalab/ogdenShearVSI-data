function out=matrix_satrecfit(varargin)
%out=matrix_satrecfit(indata,taxis[,options?])
%
%options:    
% 'fix'     :   force longest time to equal t=infinty

indata=varargin{1};
taxis=varargin{2};
taxis=(taxis(:))';

if any(strcmp('blur', varargin));
    bindex=find(strcmp('blur', varargin));
    if nargin>bindex && isnumeric(varargin{bindex+1});
        blur_range=varargin{bindex+1};
    else
        blur_range=2;
    end
    si=size(indata);
    numpics=numel(indata)/(si(1)*si(2));
    indata=reshape(indata,[si(1) si(2) numpics]);
    for jj=1:numpics;
        indata(:,:,jj)=blur(indata(:,:,jj),blur_range);
    end
    indata=reshape(indata,si);
end


if any(strcmp('fix',varargin));
	free_amp_flag=false;
else
	free_amp_flag=true;
end

if any(strcmp('fix2',varargin));
    free_amp_flag=false;
    %do a saturation revovery, and give extra weight to the longest TR data
    %point
    [maxtime,maxind]=max(taxis);
    si=size(indata);
    if ndims(indata)==3;
        pad=indata(:,:,maxind);
    else
        pad=indata(:,:,:,maxind);
    end
    %for nw=1:3;
    %    taxis=[taxis (nw/3+2)*taxis(maxind)];
    %    indata=cat(ndims(indata),indata,pad);
    %end
end


if ndims(indata) == 3;                          %one slice, multiple times
    if free_amp_flag;
        out = matrix_satrecfit_3d(indata,taxis);            %option no 'fix', 'free search for satamp
    else
        out=matrix_satrecfit_3d_t_infinity(indata,taxis);   %option 'fix'
    end
    
elseif ndims(indata) ==4;                       %multislice, multiple times
    [ntd nph nsl nti] = size(indata);
    out.amplitude=zeros(ntd,nph,nsl);
    out.Tau=zeros(ntd,nph,nsl);
    for ns=1:nsl;
        if free_amp_flag;
            oout= matrix_satrecfit_3d(squeeze(indata(:,:,ns,:)),taxis);
        else
            oout = matrix_satrecfit_3d_t_infinity(squeeze(indata(:,:,ns,:)),taxis);
        end
        out.amplitude(:,:,ns)=oout.amplitude(:,:);
        out.Tau(:,:,ns)=oout.Tau(:,:);
    end
elseif ndims(indata) == 2;
    %1. find which dimension corresponds to the time axis
    si=size(indata);
    if si(2)~= numel(taxis);
        indata=abs(indata)';
    end
    
    %create a fake 3d data set
    indata=reshape(indata,[1 size(indata)]);
    
    if free_amp_flag;
        out = matrix_satrecfit_3d(indata,taxis);
    else
        out=matrix_satrecfit_3d_t_infinity(indata,taxis);
    end
    
end

out.Tau(abs(out.Tau(:))>15)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saturation recovery fit letting amplitude and tau vary:
function out=matrix_satrecfit_3d(varargin)
%out=matrix_invrecfit2(incube,taxis)
%fit to exponential recovery
%incube: nro X npe (x numslices) X numel(taxis) is a stack of images acquired at different times
%the times passed down in the taxis.
%out is structure with fields 'amplitude' and 'tau', result of a pixelwise
%
%the inversion recovery  ' y= a*(1-exp(bt)) '
%can't be done by straight linear regression, so this routine implements an
%iterative linear regression, where an initail amplitude matrix a is
%guessed, then updated iteratively
%log[-(y-a)] = log(a) + bt; 

incube=varargin{1};
taxis=varargin{2};

%add a t=0 point, where the amplitude is zero
taxis=cat(1,0,taxis(:));
sincube=size(incube);
incube=cat(3,zeros(sincube(1:2)),incube);

[tmax,tmaxind]=max(taxis);


if ndims(incube) == 3;
    %third dimension is the echotime dimension
    phasemat=incube(:,:,tmaxind)./abs(incube(:,:,tmaxind)); %phase of the image with the longest wait time
    phincube=zeros(size(incube));
    
    %1. rotate the phase
    for nt=1:numel(taxis);
         phincube(:,:,nt)= +incube(:,:,nt)./phasemat;
    end
    
    %2. generate a starting guess for the amplitude
    a0=real(1.05*phincube(:,:,tmaxind));
    %2b estimate a good limit for the control plot
    zlim=max(a0(:))/50*[-1 1];
    
    %3. iteratively fit and update the amplitude guess
    Y=zeros(size(phincube));
    chi2=zeros(size(a0));
    director=-ones(size(a0))/50;
    niter=25;
    %hf=figure('Name','convergence monitor','Position',[1000 500 200 200]);
    %axis square
    for jj=1:niter;
        for nt=1:numel(taxis);
            Y(:,:,nt)=a0-real(phincube(:,:,nt));
        end
        oldchi2=chi2;
        [F,Tau,chi2] = matrix_fit(Y,taxis);


        dummy=size(phincube);
        plotflag=false;
        if dummy(1)>1 && plotflag;
            imagesc(chi2);
            set(gca,'clim',[0 1e8],'zlim',[0 1e8]);
            title(num2str(jj));
            pause(0.005);
        end
        
        %update the amplitude guess%%%%%
        if jj==1;
            a0 = F;
        else
            riseindex=chi2(:)>oldchi2(:);
            director(riseindex)=-director(riseindex)/2;
            a0=a0.*(1+director);
            convergenceplotflag=false;
            if convergenceplotflag;
                subplot(1,2,1);
                plot(jj,a0(100,120),'o'); hold on
                subplot(1,2,2);
                plot(jj,chi2(100,120),'o'); hold on
                pause(0.1);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %close(hf);
end

out.amplitude=a0;
out.Tau=Tau;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saturation recovery fit letting tau vary, take the longest time amplitude
% to be the t=infinity amplitude:

function out=matrix_satrecfit_3d_t_infinity(varargin)
%out=matrix_satrecfit_3d_t_infinity(varargin)
%fit to exponential recovery, fixing the amplitude using the longest time
%value
%incube: nro X npe (x numslices) X numel(taxis) is a stack of images acquired at different times
%the times passed down in the taxis.
%out is structure with fields 'amplitude' and 'tau', result of a pixelwise
%
%the inversion recovery  ' y= a*(1-exp(bt)) '
%can't be done by straight linear regression, so this routine implements an
%iterative linear regression, where an initail amplitude matrix a is
%guessed, then updated iteratively
%log[-(y-a)] = log(a) + bt; 

incube=varargin{1};
taxis=varargin{2};


%average echoes obtained with equal TR
TRaveflag=true;
if TRaveflag;
    taxistemp=taxis;
    for jj=1:numel(taxis);
        if any(taxistemp~=0);
            ind=find(taxistemp~=0, 1,'first');
            indvec = taxistemp==taxistemp(ind);
            incube(:,:,ind)=mean(incube(:,:,indvec),3);
            taxistemp(indvec)=0;
            iind(jj)=ind;
            weight(ind)=sum(indvec); %how many data there are with this time
        end
    end
    weight=weight(iind);
    taxis=taxis(iind);
    incube=incube(:,:,iind);
end

[tmax,tmaxind]=max(taxis);
[tmin,tminind]=min(taxis);


%third dimension is the echotime dimension
phasemat=incube(:,:,tmaxind)./abs(incube(:,:,tmaxind)); %phase of the image with the shortest wait time
phincube=zeros(size(incube));

phasesensitive=false;
if phasesensitive;
    %1. rotate the phase
    for nt=1:numel(taxis);
        phincube(:,:,nt)= incube(:,:,nt)./phasemat;
    end
else
    %take the absolute value
    phincube = abs(incube);
end

%2. starting guess for the amplitude
a0=abs(phincube(:,:,tmaxind));

%3. generate a new data set containing zero time and omitting infinite time
tindex = [];
for jj=1:numel(taxis);
    if jj ~= tmaxind;
        tindex=[tindex jj];
    end
end

%saturation recovery without the longest time data
newphincube=cat(3, a0/10000.*randn(size(a0)), phincube(:,:,tindex));
newtaxis=[0 taxis(tindex)];
weight= [weight(tmaxind) weight(tindex)];

%convert saturation recovery to an exponential decay, using the longest
%time data
si=size(newphincube);
for jj=1:si(3); %turn sat recovery into an exponential decay
    newphincube(:,:,jj)=abs(a0)-abs(newphincube(:,:,jj));
end
%weight(1)=10*weight(1);


%first go around, assuming longest time is infinite time
[out.amplitude out.Tau]=matrix_fit(newphincube,newtaxis,weight);

%4. now use the Tau value obtained here to get a better estimate for the
%infinite time amplitude

normalize_to_largestTR=false;
%1. first option: a0 is based on the longest time data, rescale it to account for
%finite TR
if normalize_to_largestTR;
    normalization=abs(1-exp(-tmax./abs(out.Tau)));
    normalization(normalization>1)=1;
    normalization(normalization<0.5)=0.5;
    a0=a0./normalization;
else
    %2. second option recalculate the infinte time amplitude using all echoes
    si=size(phincube);
    sum_norm=zeros(si(1:2));
    for jj=1:numel(taxis);
        Tauref=out.Tau;
        %Tauref(Tauref(:)>4)=2;              %set large T1 outliers (in the noise, mainly) to 2s
        %Tauref(Tauref(:)<min(taxis)/10)=2;   %set unmeasurable short T1s to 2s
        normalization=abs(1-exp(-taxis(jj)./abs(blur(Tauref,1))));
        normalization(normalization>1)=1;
        normalization(normalization<0.2)=0.2;
        sum_norm=sum_norm+normalization;
    end
    new_a0=sum(phincube,3)./sum_norm;
    a0=new_a0;
end

%initialize new data
newphincube=cat(3, zeros(size(a0)), phincube(:,:,tindex));
si=size(newphincube);
for jj=1:si(3); %turn sat recovery into an exponential decay
    newphincube(:,:,jj)=abs(a0)-abs(newphincube(:,:,jj));
end
%second go around, assuming longest time is infinite time
[out.amplitude out.Tau]=matrix_fit(newphincube,newtaxis);



function [amplitude, tau, chi2] = matrix_fit(varargin)
incube=varargin{1};
maxincube=max(incube(:));
incube=incube/maxincube;
taxis=varargin{2};
if nargin==3;
    weights=varargin{3};
else
    weights=ones(size(taxis));
end

%this routine fits the log of the signal to a line
%log(y) = log(a) + bt
noiselevel=0.1;
oneovers2=zeros(size(incube));
for jj=1:numel(taxis);
    oneovers2(:,:,jj) = weights(jj)*(abs(incube(:,:,jj))/noiselevel).^2; %scale noise to account for change of variables to log(y);
end
%oneovers2(incube(:)<0)=10;
%incube(incube<0)=1e-4;

sincube=size(incube);
lincube=log(abs(incube));

x=ones(sincube);
for jj=1:sincube(3);
    x(:,:,jj)=x(:,:,jj)*taxis(jj);
end
y=lincube;
xy=x.*y;
x2=x.*x;

Sx2 =sum(x2.*oneovers2,3);
Sy  =sum(y.*oneovers2,3);
Sx  =sum(x.*oneovers2,3);
Sxy =sum(xy.*oneovers2,3);
Soneovers2=sum(oneovers2,3);

Delta=Soneovers2.*Sx2-Sx.^2;
b= (Soneovers2.*Sxy - Sx.*Sy)./Delta;
a= (Sx2.*Sy-Sx.*Sxy)./Delta;

amplitude=exp(a)*maxincube;
tau=-1./b;
R=-b;


%errors in the fitted parameters:
chi2cube=zeros(size(incube));
for jj=1:numel(taxis);
    chi2cube(:,:,jj)=((amplitude.*exp(b*taxis(jj))-maxincube*incube(:,:,jj))).^2;
end

chi2=sum(chi2cube,3);


    