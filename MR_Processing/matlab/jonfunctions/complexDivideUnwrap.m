%Uli's code


%now calculate the numeric drivatives for the strain image

%Note: complex image is the 256 x 48 x 40 image, i.e. doubled size in the
%phase encode (PE) and slice (SL) dimensions
CI = cat(4,hCIx,hCIy,hCIz);

si=size(CI);
complex_strainimage=zeros([si 3]);
CSI = complex_strainimage;
ind1=[2:si(1)-1];
ind2=[2:si(2)-1];
ind3=[2:si(3)-1];

phaseimage=CI./abs(CI);

% RO encode            lambda  subtraction
complex_strainimage(ind1,:,:,1,1)=(phaseimage(ind1+1,:,:,1)./phaseimage(ind1-1,:,:,1));
complex_strainimage(:,ind2,:,1,2)=(phaseimage(:,ind2+1,:,1)./phaseimage(:,ind2-1,:,1));
complex_strainimage(:,:,ind3,1,3)=(phaseimage(:,:,ind3+1,1)./phaseimage(:,:,ind3-1,1));

% PE encode          lambda  subtraction
complex_strainimage(ind1,:,:,2,1)=(phaseimage(ind1+1,:,:,2)./phaseimage(ind1-1,:,:,2));
complex_strainimage(:,ind2,:,2,2)=(phaseimage(:,ind2+1,:,2)./phaseimage(:,ind2-1,:,2));
complex_strainimage(:,:,ind3,2,3)=(phaseimage(:,:,ind3+1,2)./phaseimage(:,:,ind3-1,2));

% SL strain         lambda  subtraction
complex_strainimage(ind1,:,:,3,1)=(phaseimage(ind1+1,:,:,3)./phaseimage(ind1-1,:,:,3));
complex_strainimage(:,ind2,:,3,2)=(phaseimage(:,ind2+1,:,3)./phaseimage(:,ind2-1,:,3));
complex_strainimage(:,:,ind3,3,3)=(phaseimage(:,:,ind3+1,3)./phaseimage(:,:,ind3-1,3));

%%
%Jon's version using a weighted matrix
%Note: have to correct the weighting for the linear coordinate being multiplied in
x = [2 1 0 1 2];
d = [0.109603762960256,0.276690988455550,0,-0.276690988455550,-0.109603762960256]./x; d(isnan(d)) = 0;
d = d/(sum(d(1:floor(length(d)/2)).*x(1:floor(length(d)/2))));
wlen = length(d); remsz = (wlen-1)/2;
ind1=(1+remsz):(si(1)-remsz);
ind2=(1+remsz):(si(2)-remsz);
ind3=(1+remsz):(si(3)-remsz);
for i=1:3
    CSI(ind1,:,:,i,1)=((phaseimage(ind1+1,:,:,i)./phaseimage(ind1-1,:,:,i))).^d(2).*((phaseimage(ind1+2,:,:,i)./phaseimage(ind1-2,:,:,i))).^d(1);
    CSI(:,ind2,:,i,2)=((phaseimage(:,ind2+1,:,i)./phaseimage(:,ind2-1,:,i))).^d(2).*((phaseimage(:,ind2+2,:,i)./phaseimage(:,ind2-2,:,i))).^d(1);
    CSI(:,:,ind3,i,3)=((phaseimage(:,:,ind3+1,i)./phaseimage(:,:,ind3-1,i))).^d(2).*((phaseimage(:,:,ind3+2,i)./phaseimage(:,:,ind3-2,i))).^d(1);
end


%Voxel lengths delta
try
    delta(1)=HIRES.axis1(2)-HIRES.axis1(1);
    delta(2)=HIRES.axis2(2)-HIRES.axis2(1);
    delta(3)=HIRES.axis3(2)-HIRES.axis3(1);
catch
    delta(1)=axis1(2)-axis1(1);
    delta(2)=axis2(2)-axis2(1);
    delta(3)=axis3(2)-axis3(1);
end

%%
GPhij = zeros([si 3]);
tic
for m=1:3
    tic
    [GPhij(:,:,:,m,1), GPhij(:,:,:,m,2), GPhij(:,:,:,m,3)] = complexDivideN(phaseimage(:,:,:,m),'optimal5');
    toc
end
toc

strainimage=angle(complex_strainimage);
SSI = angle(CSI);
GUij = angle(GPhij);
for e =1:3
    for g =1:3
        strainimage(:,:,:,e,g)=...
            strainimage(:,:,:,e,g)/(2*pi)*lambda(e)/(2*delta(g));
        SSI(:,:,:,e,g)=...
            SSI(:,:,:,e,g)/(2*pi)*lambda(e)/(2*delta(g));
        GUij(:,:,:,e,g)=...
            GUij(:,:,:,e,g)/(2*pi)*lambda(e)/(2*delta(g));
    end
end

newmask1 = mask;
%histogram(strainimage(125:175,:,:,1,1).*newmask1(125:175,:,:)); hold on;
view3dgui(-GUij(:,:,:,1,1).*newmask1);
figure; 
%z=60;
for p=1:3
    for q = 1:3
        Gjui{p,q} = -GUij(:,:,:,p,q);
        figure(301);
        histogram(-GUij(125:175,:,:,p,q).*newmask1(125:175,:,:)); hold on;
        figure; 
        subplot(1,3,2)
        imagesc(-GUij(:,:,z,p,q).*newmask1(:,:,z)); axis image; set(gca,'ydir','normal'); colormap magma;
        subplot(1,3,1)
        imagesc(strainimage(:,:,z,p,q).*newmask1(:,:,z)); axis image; set(gca,'ydir','normal'); colormap magma;
        subplot(1,3,3)
        imagesc((-GUij(:,:,z,p,q)-strainimage(:,:,z,p,q)).*newmask1(:,:,z)); axis image; set(gca,'ydir','normal'); colormap magma;        
    end
end
histogram(SSI(125:175,:,:,1,2).*newmask1(125:175,:,:))

epsij = calculateEpsij(Gjui);
Eij = calculateEij(Gjui);
for p=1:3
    for q = 1:3
        figure(302);
        histogram(Eij{p,q}(125:175,:,:).*newmask1(125:175,:,:)); hold on;
        figure(401); 
        subplot(3,3,3*(q-1)+p)
        imagesc(Eij{p,q}(:,:,z).*newmask1(:,:,z)); axis image; set(gca,'ydir','normal'); colormap magma;
    end
end

% figure;
% subplot(1,3,1);
% imagesc(Eij{1,1}(:,:,60).*mask(:,:,60)); axis image; set(gca,'ydir','normal'); colormap magma;
% subplot(1,3,2);
% imagesc(Eijnew{1,1}(:,:,60).*mask(:,:,60)); axis image; set(gca,'ydir','normal'); colormap magma;
% subplot(1,3,3);
% imagesc((Eij{1,1}(:,:,60)-Eijnew{1,1}(:,:,60)).*mask(:,:,60)); axis image; set(gca,'ydir','normal'); colormap magma;


Em = calculateEm(Eij);
figure(402);
imagesc(Em(:,:,z).*newmask1(:,:,z)); axis image; set(gca,'ydir','normal'); colormap magma;
%%
%newmask1(log10(abs(CI(:,:,:,1)))<1.5) = nan;
% figure; view3dgui(abs(CI(:,:,:,1)))
% hold on; view3dgui(SSI(:,:,:,1,1).*newmask1)

% 
% DESTE.strains=strainimage;
% DESTE.lambda=lambda/1000;   %in mm
% DESTE.comment=wigglecomment;
% DESTE.timest=timestampvec;
% 
% 
% 
% %%
% 
% 
% p = repmat(1:256,[48 1])';
% gp = convn(p,[-1; 0; 1],'same')/2;
% 
% phix = [zeros(1,16) linspace(0,10*pi,256-32) ones(1,16)*10*pi]; phix = repmat(phix,[48 1])';
% %phiy = linspace(0,0,48); phiy = repmat(phiy,[256 1]);
% zx = cos(phix) + 1i*sin(phix);
% ax = real(zx); bx = imag(zx);
% ang_zx = angle(zx);
% imagesc(ang_zx);
% conv_ang_zx = convn(ang_zx,[-1; 0; 1],'same');
% imagesc(conv_ang_zx);



