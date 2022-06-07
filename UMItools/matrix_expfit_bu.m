function out=matrix_expfit(varargin)
%fits an exponential decay to a 3d or 4d data set, last dimension is time
%out=matrix_expfit(incube,taxis [,'oddeven'],'weighting', weight matrix)
%incube: size(nro npe numel(taxis)) is a stack of single slice images acquired at different times
%
%the option oddeven takes the running average 
%out is structure with fields 'amplitude' and 'tau', 'chi2', result of a
%pixelwise fit

incube=varargin{1};
taxis=varargin{2};
[ny nx ntimes]=size(incube);
if ntimes~=numel(taxis);
    display('matrix_expfit: timeaxis has to match number of echoes in incube!')
    return
end

%sort by ascending times
[taxis, indvec]=sort(taxis);
if ndims(incube)==3;
    incube=incube(:,:,indvec);
elseif ndims(incube)==4;
    incube=incube(:,:,:,indvec);
end


if any(strcmp(varargin,'weighting'));
    dummy=find(strcmp(varargin,'weighting'));
    if nargin>dummy;
        if isequal(size(varargin{dummy+1}),size(incube));
            weightcube=varargin{dummy+1};
            if ndims(incube)==3;
                weightcube=weightcube(:,:,indvec);
            elseif ndims(incube)==4;
                weightcube=weightcube(:,:,:,indvec);
            end
       else
           display('weighting in expfit must have same size as data');
       end
   else
       display('weighting indicated to matrix_expfit, but no weighting passed down');
   end
end

%option to fit to an average of adjacent echoes, to smooth out oddeven variability %
if any(strcmp(varargin, 'oddeven'));
    %diagnostics whether off even asymmetry is an issue
    diagflag=false;
    if diagflag;
        si=size(incube);
        etl=si(end);
        testm=sum(incube,ndims(incube)); %sums over the last dimension
        if ndims(testm)==2; %matrix, single slice multiecho
            nsli=1;
            mask=false(si(1:2));
        else
            si=size(testm);
            nsli=si(3);
            mask=false(si(1:3));
        end
        
        difference=zeros(etl,nsli);
  
            
        for sli=1:nsli;
            mask=testm(:,:,sli)>0.5*max(max(testm(:,:,sli))); %mask out everything that is halfthe max amplitude
            %go over all echoes, get a mean amplitude over the mask  ******
            for ne=1:etl;
                if nsli==1;
                    patch=incube(:,:,ne);
                else
                    patch=incube(:,:,sli,ne);
                end
                amp(ne,sli)=sum(patch(mask(:)));
            end
            % *************************************************************
            
            for ne=2:etl-1;
                difference(ne,sli)=abs((amp(ne-1,sli)+amp(ne+1,sli))/2-amp(ne,sli))/((amp(ne-1,sli)+amp(ne+1,sli))/4+amp(ne,sli)/2);
            end
            
        end
        
    end
    
    newincube=zeros(ny,nx,ntimes-1);
    newtaxis=zeros(1,ntimes-1);
    for jj=1:ntimes-1;
        newincube(:,:,jj)=(incube(:,:,jj)+incube(:,:,jj+1))/2;
        newtaxis(jj)=(taxis(jj)+taxis(jj+1))/2;
    end
    incube=newincube;
    taxis=newtaxis;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%noiselevel
if exist('weightcube', 'var');
    noiselevel=abs(weightcube+10^6*eps);
else
    noiselevel=estimate_noiselevel(incube)*ones(size(incube)); 
end

%single slice multiecho
if ndims(incube) == 3;
    %third dimension is the echotime dimension
    [amplitude, Tau, chi2]=matrix_fit(abs(incube),taxis,noiselevel);
    nosl=1;
end

%multislice multiecho
if ndims(incube) == 4;
    %third dimension is the slice dimension
    %fourth dimension is the echo time dimension
    sincube=size(incube);
    nosl=sincube(3);
    amplitude=zeros(sincube(1:3));
    Rate=(zeros(sincube(1:3)));
    for ns=1:nosl;
        [amplitude(:,:,ns), Tau(:,:,ns), chi2(:,:,ns)] = ...
            matrix_fit(squeeze(incube(:,:,ns,:)),taxis,squeeze(noiselevel(:,:,ns,:)));
    end
end

% option to mask out noisy data (cosmetic) ********************************
find95th=false;
if find95th;
    %find the 95th percentile amplitude;
    MA=max(amplitude(:));
    normamp=abs(amplitude(:)/MA);
    limitvec=0:0.05:1;
    N=histc(normamp,limitvec);
    N(1)=0; %null the noise floor
    N(2)=0; %null the noise floor
    nof_finite_pixels=sum(N); %how many finite pixels are there?
    integral=cumsum(N)/nof_finite_pixels;
    percentile95_index=find(integral>0.95,1,'first');
    level95=limitvec(percentile95_index)*MA;
    
    m%ask as noise all pixels whose level is than 5% of the 95% level
    noiseindex=abs(amplitude(:))<0.05*level95;
    amplitude(noiseindex)=0;
    Tau(noiseindex)=0;
end
%**************************************************************************

out.amplitude=amplitude;
out.Tau=Tau;
out.chi2=chi2;


function [amplitude, tau, chi2] = matrix_fit(incube,taxis,noiselevel)
%this routine fits the log of the signal to a line
%log(y) = log(a) + bt
oneovers2 = (abs(incube)./noiselevel).^2; %scale noise to account for change of variables to log(y);
%oneovers2 = ones(size(oneovers2));
sincube=size(incube);
y=log(abs(incube));

x=ones(sincube);
for jj=1:sincube(3);
    x(:,:,jj)=x(:,:,jj)*taxis(jj);
end

xy=x.*y;
x2=x.*x;

Sx2=sum(x2.*oneovers2,3);
Sy=sum(y.*oneovers2,3);
Sx=sum(x.*oneovers2,3);
Sxy=sum(xy.*oneovers2,3);
Soneovers2=sum(oneovers2,3);

Delta=Soneovers2.*Sx2-Sx.^2;
b= (Soneovers2.*Sxy - Sx.*Sy)./Delta;
a= (Sx2.*Sy-Sx.*Sxy)./Delta;

amplitude=exp(a);
tau=-1./b;

%chi2 pixel by pixel
synthdata=zeros(size(incube));
for jj=1:sincube(3);
    synthdata(:,:,jj)=amplitude.*exp(b.*x(:,:,jj));
end
chi2=sum((incube-synthdata).^2,3);

function noiselevel=estimate_noiselevel(incube)
%estimate the noise using NxN patches of the last echo
N=10;
lastecho=incube(:,:,1);
si=size(lastecho);
ni=floor(si(1)/N);
nj=floor(si(2)/N);
noisetest=zeros(ni,nj);
for ii=0:ni-1;
    for jj=0:nj-1;
        ivec=(1:N) + ii*N;
        jvec=(1:N) + jj*N;
        patch=lastecho(ivec,jvec);
        noisetest(ii+1,jj+1)=sqrt(var(real(patch(:)))+var(imag(patch(:))));
    end
end
noiselevel=mean(noisetest(:));




