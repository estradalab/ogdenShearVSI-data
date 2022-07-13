function out=matrix_invrecfit_oddeven(indata,taxis,varargin)
%out=matrix_invrecfit_up_to_zerocrossing(incube,taxis,weighting)
%incube is a stack of image matrices (the same slice) acquired at different times
%taxis is an axis of times;

options=varargin;

if ndims(indata) == 3;                          %one slice, multiple times
    out = matrix_invrecfit_3d(indata,taxis,options);
elseif ndims(indata) ==4;                       %multislice, multiple times
    [ntd nph nsl nti] = size(indata);
    out.amplitude=zeros(ntd,nph,nsl);
    out.Tau=zeros(ntd,nph,nsl);
    for ns=1:nsl;
        oout= matrix_invrecfit_3d(squeeze(indata(:,:,ns,:)),taxis,options);
        out.amplitude(:,:,ns)=oout.amplitude(:,:);
        out.T1(:,:,ns)=oout.Tau(:,:);
    end
end
out.T1(abs(out.T1(:))>5)=0;

function out=matrix_invrecfit_3d(incube,taxis,varargin)
%out=matrix_invrecfit_3d(incube,taxis,weighting)
%fit to exponential recovery
%incube: nro X npe (x numslices) X numel(taxis) is a stack of images acquired at different times
%the times passed down in the taxis.
%out is structure with fields 'amplitude' and 'T1', result of a pixelwise
%
%the inversion recovery  ' y= a*(1-2*exp(bt)) '
%can't be done by straight linear regression, so this routine implements an
%iterative linear regression, where an initial amplitude matrix a is
%guessed, then updated iteratively

options=varargin;

[taxis,sortindex]=sort(taxis); %sort by ascending times
incube=incube(:,:,sortindex);  %sort by ascending times

[tmax,tmaxind]=max(taxis);
[tmin,tminind]=min(taxis);

%Third dimension is the echotime dimension
phase_matrix=incube(:,:,tminind)./abs(incube(:,:,tminind)); %phase of the image with the shortest wait time
phased_incube=zeros(size(incube));

%1. rotate the phase, short times are positive, long times negative
for nt=1:numel(taxis);
    phased_incube(:,:,nt)= incube(:,:,nt)./phase_matrix; %rotate signal into the imaginary axis
end

%2. generate a starting guess for the amplitude, using the two shortest
%times and interpolating backwards
a0= abs(phased_incube(:,:,1) - (phased_incube(:,:,2)-phased_incube(:,:,1))/(taxis(2)-taxis(1)) * taxis(1));

%3. iteratively fit and update the amplitude guess
si=size(phased_incube);
Y=zeros(size(phased_incube));

niter=20;
for jj=1:niter;
    for nt=1:numel(taxis);
        Y(:,:,nt)= phased_incube(:,:,nt)+a0;
    end
    outY=matrix_expfit(Y,taxis,options);
    
    F=outY.amplitude/2;
    mesh(F-a0);
    set(gca,'clim',5e3*[-1 1],'zlim',[-5000 5000]); title(num2str(jj)); pause(0.1);
    %imagesc(gg);
    
    %update the amplitude guess%%%%%
    if jj==1;
        a0 = F/2;
    else
        a0=(a0+F)/2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

out.amplitude=F;
out.Tau=outY.Tau;


