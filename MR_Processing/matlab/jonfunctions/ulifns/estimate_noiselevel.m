function [signallevel, noiselevel]=estimate_noiselevel(varargin)
incube=varargin{1};
%estimate the noise using NxN patches in dimensions 1 & 2 of the dataset
%default 8x8 pixel patch
%default estimate for empty space in the data set = 25%
N=8;                    %the size of a patch of pixels used to coarse-grain the image;
emptyspacefraction=0.25;%estimate of how much empty space there is in the image data, for the noise calculation

si=size(incube);
incube=reshape(incube,[si(1) si(2) numel(incube)/(si(1)*si(2))]);
refplane=incube(:,:,(1 : end));
si=size(refplane);
ni=floor(si(1)/N);
nj=floor(si(2)/N);
nk=floor(si(3)/N);
noisetest=zeros(ni,nj,nk);
amplitudetest=zeros(ni,nj,nk);
for kk=0:nk-1;
    for ii=0:ni-1;
        for jj=0:nj-1;
            ivec=(1:N) + ii*N;
            jvec=(1:N) + jj*N;
            kvec=(1:N) + kk*N;
            patch=refplane(ivec,jvec,kvec);
            noisetest(ii+1,jj+1,kk+1)=sqrt(var(abs(patch(:))));
            amplitudetest(ii+1,jj+1,kk+1)=abs(mean(patch(:)));
        end
    end
end

sortednoiselevel=sort(noisetest(:));
sortedamplitudelevel=sort(amplitudetest(:));

% the noise level is calculated where the signal is zero,
% (or finite for abs value data
% by taking the mean of the first 25% percent of noise (asumed empty)
% and the mean of the first 25% of the signal
noiselevel= mean(sortednoiselevel(1:round(numel(sortednoiselevel)*emptyspacefraction))) + ...
    mean(sortedamplitudelevel(1:round(numel(sortedamplitudelevel)*emptyspacefraction)));


%signal level is calculate at the end of the signal distribution
signallevel=mean(sortedamplitudelevel(round(end*0.95):(end-10)));

if any(strcmp(varargin,'plot'));
    figure;
    plot(sortedamplitudelevel,'r');
    hold on;
    plot(sortednoiselevel,'b');
    set(gca,'Yscale','log');
    legend({'signal'; 'noise'});
    legend boxoff;
    title(['S=' num2str(round(signallevel)) ' / N=' num2str(round(noiselevel))]);
end
