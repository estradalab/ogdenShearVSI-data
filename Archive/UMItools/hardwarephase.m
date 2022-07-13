function [out, meanSEphi]=hardwarephase(varargin)
%[out,slicephase]=hardwarephase(in, [echoindex])
%in is an image structure returned by varianms or an image matrix
%out is an image structure or image matrix, phase corrected with a common
%phase plane to take care p0 and p1 in the readout direction.
%if there are multiple echoes, the echo index points to the true spin echo,
%where field inhomogeneities are perfectly refocussed.
%This reference echo will only have the arbitrary common phase p0, and a
%linear phase in the readout direction (1st dimension) due to off-center
%echoes.


in=varargin{1};
if nargin>1;
	echoindex=varargin{2};
else
    display('Spin echo index defaulted to 2nd echo in hardwarephase.m ...');
    echoindex=2;
end

instructflag=isstruct(in);

if instructflag;
	out=in;
	si=size(in.image);
	if ndims(in.image)==4; %multislice multi-echo
		nofslices=si(3);
		for ns=1:nofslices;
			spinecho=in.image(:,:,ns,echoindex);       %this is the echo without phase, except the hardware phase
			biggestecho=in.image(:,:,ns,1);            %this is the first echo, used to construct a mask
			[phasepl,meanSEphi(ns)]=phaseplane(spinecho,biggestecho);
			for ne=1:si(4);
				out.image(:,:,ns,ne)=in.image(:,:,ns,ne).*exp(-1i*phasepl);
			end
        end
        
        %all the images have been phase adjusted to the offset phases of
        %the spin echo; they could still have a phase due to offset
        %frequency and timing
        %for ns=1:nofslices;
        %    for ne=1:si(4);
        %        [dummy,meanPH(ns,ne)]=phaseplane(out.image(:,:,ns,ne),biggestecho);
        %    end
        %end
        
	elseif ndims(in.image)==3; %single slice multi-echo
        nofslices=1;
		spinecho=in.image(:,:,echoindex);       %this is the echo without phase, except the harware phase
		biggestecho=in.image(:,:,1);            %this is the first echo, used to construct a mask
		[phasepl,meanSEphi(1)]=phaseplane(spinecho,biggestecho);
		for ne=1:si(4);
			out.image(:,:,ne)=in.image(:,:,ne).*exp(-1i*phasepl);
		end
    end
else
    %hardware phase on a single echo image
	spinecho=in;
	biggestecho=in;
	[phasepl,meanSEphi(1)]=phaseplane(spinecho,biggestecho);
	out=spinecho*exp(-1i*phasepl);
end

function [phasepl, meanSEphase]=phaseplane(spinecho,biggestecho)
%1. make the mask
si=size(spinecho);
Msk=false(size(biggestecho));
magn=abs(blur(abs(biggestecho),2));
magn=magn/max(magn(:));
limvec=(0:0.05:1);
ct=histc(magn(:),limvec);
ct(1)=0; %reject empty space and 5% above it
ct(2)=0; %reject 5% to 10%

%Now find those parts of the distribution where a threshold percentage
%of the signal resides;
thr=0.95;
cumct=cumsum(ct)/sum(ct);                 %cumulative sum
maxind=find(cumct<thr,1,'last');
Msk=magn>0.3*limvec(maxind);              %mask for each slice, based on strongest echo (i.e first)

anglespinecho=angle(spinecho);

%get the y slope
diffy=diff(anglespinecho,1,1);
diffy=cat(1,diffy,diffy(end,:));
%mask out the jumps
diffy(abs(diffy)>pi)=NaN;
diffy(~Msk(:))=NaN;
[fittedplane,gx,gy]=planefit(diffy);
ymatrix=gy*[1:si(1)]'*ones(1,si(2));
ymatrix=ymatrix-mean(ymatrix(:));

%correct for the slope
spinecho_detilt=spinecho.*exp(-1i*ymatrix);
spinecho_detilt(~Msk(:))=NaN;

%get the mean phase of the image, account for wrap around
stdph(1)=std(angle(spinecho_detilt(~isnan(spinecho_detilt))));
stdph(2)=std(angle(spinecho_detilt(~isnan(spinecho_detilt))*exp(1i*pi/2)));
stdph(3)=std(angle(spinecho_detilt(~isnan(spinecho_detilt))*exp(1i*pi)));
stdph(4)=std(angle(spinecho_detilt(~isnan(spinecho_detilt))*exp(1i*3*pi/2)));

mnph(1)=mean(angle(spinecho_detilt(~isnan(spinecho_detilt))));
mnph(2)=mean(angle(spinecho_detilt(~isnan(spinecho_detilt))*exp(1i*pi/2)))-pi/2;
mnph(3)=mean(angle(spinecho_detilt(~isnan(spinecho_detilt))*exp(1i*pi)))-pi;
mnph(4)=mean(angle(spinecho_detilt(~isnan(spinecho_detilt))*exp(1i*3*pi/2)))-3*pi/2;

mnph=mod(mnph,2*pi);

[dummy,loweststdindex]=min(stdph);
phasepl=ymatrix+mnph(loweststdindex);
meanSEphase=mnph(loweststdindex);

meanSEphase(meanSEphase>pi)=meanSEphase(meanSEphase>pi)-2*pi;


