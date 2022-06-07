function out=lig_test_averaging(P,blurvec,varargin)
% calculate the diagonal strains by shift and subtract
% subtract the mean strain, Fourier transform, and get the amplitude and
% frequency of wiggles in each direction

%% this is where the filtering and upsampling occurs
gridvec=[192 64 64];
compleximage=zeros([gridvec 3]);

orig_compleximage=P.origdata;  %this is unblurred (!) stimulated echo data, combined + and minus encoding


for nf=1:3;
     %create a temporarystucture with fields image and parss, to blur using mm voxel sizes
	temp.image=orig_compleximage(:,:,:,nf);
	temp.pars=P.origpars(1);
	temp.acquisitiontype='2D';
	
	% blur the STE with the required blurring, before upsampling
	temp=blur3d(temp,'vox', blurvec);
	display(['blur size= ' num2str(blurvec) 'mm']); %#ok<*DISPLAYPROG>
    
	% upsample without blurring
    temp=blur3d(temp,'vox',[0 0 0],'grid',gridvec);
	
	%update the complex image and phase factor to the new grid size
	compleximage(:,:,:,nf)=temp.image;
end




%% now calculate the numeric drivatives for the strain image
si=size(compleximage);
complex_strainimage=zeros([si 3]);
ind1=[2:si(1)-1];
ind2=[2:si(2)-1];
ind3=[2:si(3)-1];

phaseimage=compleximage./abs(compleximage);

% RO encode            lambda  subtraction
complex_strainimage(ind1,:,:,1,1)=(phaseimage(ind1+1,:,:,1)./phaseimage(ind1-1,:,:,1));
complex_strainimage(:,ind2,:,1,2)=(phaseimage(:,ind2+1,:,1)./phaseimage(:,ind2-1,:,1));
complex_strainimage(:,:,ind3,1,3)=(phaseimage(:,:,ind3+1,1)./phaseimage(:,:,ind3-1,1));
% PE encode          lambda  subtraction
complex_strainimage(ind1,:,:,2,1)=(phaseimage(ind1+1,:,:,2)./phaseimage(ind1-1,:,:,2));
complex_strainimage(:,ind2,:,2,2)=(phaseimage(:,ind2+1,:,2)./phaseimage(:,ind2-1,:,2));
complex_strainimage(:,:,ind3,2,3)=(phaseimage(:,:,ind3+1,2)./phaseimage(:,:,ind3-1,2));
% SL strain         lambda  subtraction
complex_strainimage(:,:,ind3,3)=(phaseimage(:,:,ind3+1,3)./phaseimage(:,:,ind3-1,3));

complex_strainimage(ind1,:,:,3,1)=(phaseimage(ind1+1,:,:,3)./phaseimage(ind1-1,:,:,3));
complex_strainimage(:,ind2,:,3,2)=(phaseimage(:,ind2+1,:,3)./phaseimage(:,ind2-1,:,3));
complex_strainimage(:,:,ind3,3,3)=(phaseimage(:,:,ind3+1,3)./phaseimage(:,:,ind3-1,3));

delta(1)=P.axis1(2)-P.axis1(1);
delta(2)=P.axis2(2)-P.axis2(1);
delta(3)=P.axis3(2)-P.axis3(1);

strainimage=angle(complex_strainimage);
lambda=P.lambda;
for encoding_directions =1:3;
	for gradient_directions =1:3;
		strainimage(:,:,:,encoding_directions,gradient_directions)=...
			strainimage(:,:,:,encoding_directions,gradient_directions)/(2*pi)*lambda(encoding_directions)...
            /(2*delta(gradient_directions)).*P.mask_hr;
	end
end

%% treat the strain image in a slice of interest, so that the fft is maximally meaningful for 

sli=32;
if any(strcmp(varargin,'encdir'));
    ind=find(strcmp(varargin,'encdir'));
    encdir=varargin{ind+1};
else
    encdir=1;
end


testslice=strainimage(:,:,sli,encdir,encdir);
testmask=P.mask_hr(:,:,sli);

insample=testmask~=0;
meanstrain=mean(testslice(insample));

%testslice(insample)=testslice(insample)-meanstrain;

%eliminate the initial couple of pixels
testprofile=sum(testmask,2);
testprofile=testprofile/max(testprofile);
insamp1=find(testprofile>0.5);
center1=mean(insamp1);
width1=abs(insamp1(1)-insamp1(end));

dummy1=1:numel(testprofile);
dummy1=dummy1-center1;
profile1=exp(- dummy1.^6/(2*(width1/2.5)^6)); profile1=profile1/max(profile1);
si=size(testslice);
profile2d=profile1(:)*ones(1,si(2));



testslice=testslice.*profile2d;

dummy=fft(testslice,2*si(1));
kmax=1/(P.origpars(1).lro*10/si(1));
kaxis=(1:2*si(1))/(2*si(1))*kmax;

spectrum=sum(abs(dummy(1:si(1),28:36)),2);

if ~any(strcmp(varargin,'noplot'))
    if any(strcmp(varargin,'overlay'))
        hold on;
    else
        figure;
    end
    plot(kaxis(1:si(1)),log10(spectrum));
    %set(gca,'Ylim',[-2 6]);
    grid on;
    title(['blurvec=' num2str(blurvec(1))  ',' num2str(blurvec(2)) ',' num2str(blurvec(3)), ';  \lambda=' num2str(lambda(1)) ',' num2str(lambda(2)) ',' num2str(lambda(3))]);
    
end







out.strainimage=strainimage;
out.blurvec=blurvec;
out.lambda=P.lambda;
out.axis1=P.axis1;
out.axis2=P.axis2;
out.axis3=P.axis3;
out.kaxis=kaxis(1:si(1));
out.spectrum=spectrum;
