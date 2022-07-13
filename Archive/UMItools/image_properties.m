function out=image_properties(varargin)
%out=image_properties(inimage [,options]);

M=varargin{1};

ca=strcmp('cut',varargin) | strcmp('cutoff',varargin);
if any(ca);
    cutoff = varargin{find(ca)+1};
    M(M>cutoff)=0;
end

%1. get some basic statistics about the image
maxM=max(M(:));                 %Maximum
bins=maxM*(0:0.05:1);           %bins for a histogram
histM=histc(M(:),bins);         %histogram
meanM=mean(M(M(:)>0.05*maxM));  %mean of finite values (>5%)

%2. make an amplitude mask
Mask=M>meanM/10;                %blank the lowest 

%high pass filter,'cut',2.5

Mprime=M;
Mprime(~Mask)=meanM;
bfvec=0.8:1:20;
for bf=1:numel(bfvec);
    B1=blur(Mprime,bfvec(bf));
    difN=Mprime-B1;
    difN(~Mask)=NaN;
    
    numberofvalidpix=sum(~isnan(difN(:)));
    valid_difN=difN(Mask(:));
    valid_absdifN=abs(valid_difN);
    [orderedNoise,oNvec]=sort(valid_absdifN);
    lowsigma(bf)=orderedNoise(floor(2*numberofvalidpix/3)); %assuming gaussian noise, 60% within pm sigma
    sigma(bf)=std(difN(Mask(:)));
    maxDeviation(bf)=orderedNoise(end);
end

out.mean=meanM;
out.SNR=sigma/meanM;
out.bottomSNR=lowsigma/meanM;
out.maxDeviation=maxDeviation/meanM;
