function output=deste_fid(timestamp,varargin)
% displacement-encoding-stimulated echo
% out=deste(timestamp,varargin)
% script to read the amplitude and phase of stimulated echo images
% for the ligament/tendon project

options=[varargin 'noiso']; %suppress the isotropic forcing

%1. identify the files with the time stamp
if isnumeric(timestamp);
    timestamp=num2str(timestamp);
end

if numel(timestamp)~=4;
    display('Identify your data using the 4 digit time stamp!') 
    return
end

d=dir;
time_flag=false(numel(d),1);
img_flag=time_flag;
PH_flag=time_flag;
fid_flag=time_flag;

for jj=1:numel(dir);
    filename=d(jj).name;
    [fpath,fname,fext]=fileparts(filename);
    
    %get or fake the last three characters of the filename to identfy phase
    %files 
    if numel(fname)>2;
        fnamelast3characters=fname((length(fname)-2):length(fname));
    else
        fnamelast3characters='xxx';
    end
    
    time_flag(jj)=~isempty(strfind(fname,timestamp));
    img_flag(jj)=~isempty(strfind(fext,'img'));
    PH_flag(jj)=strcmp(fnamelast3characters,'_PH');
    fid_flag(jj)=~isempty(strfind(fext,'fid'));
end

fid_directory=d(time_flag & fid_flag).name;
img_directory=d(time_flag & img_flag).name;
phi_directory=d(time_flag & PH_flag).name;


%get and reconstruct the multislice image from the fid directory
%the image is 4-d three spatial dimensions and plus minus lambda.

out=varianms(fid_directory, 'noiso','flipfb'); %straight up reconstruction from fid

phasefactor=out.image(:,:,:,1)./out.image(:,:,:,2);
phasefactor=phasefactor./abs(phasefactor);

compleximage=mean(abs(out.image),4).*phasefactor;

blurflag=true;
if blurflag;
    if any(strcmp(varargin,'grid'));
        ind=find(strcmp(varargin,'grid'));
        gridvec=varargin{ind+1};
    else
        gridvec=size(compleximage);
    end
    %blur 3d
    dummy.image=compleximage;
    dummy.pars=out.pars;
    dummy=lig_blur3d(dummy,'vox',0.6*[1 1 1],'grid',gridvec);
    compleximage=dummy.image;
end


    

%get the paramters out of the fid directory
pars=getparameters([fid_directory '/procpar']);

%construct the complex image
output.image=compleximage;

%extract the required parameters for plotting, construct the axes
%all lengths in cm
[numro numpe numsl]=size(output.image);
output.pars.xaxis=[0:numpe-1]/numpe*pars.lpe;
output.pars.yaxis=[0:numro-1]/numro*pars.lro;
output.pars.zaxis=pars.pss;
output.pars.lambda=pars.lambda;

%claculate the diffusive attenuation factor b=k^2 Delta
%standard units are s/mm^2
output.pars.diffusive_b=(2*pi/(pars.lambda*10))^2*pars.TS;


if any(strcmp(options,'plot'));
    pic=compleximage;
    
    %plot the slices
    si = size(pic);
    figure('position' , [50 50 1200 400]);
    ns=si(3);
    
    for sl=1:ns;
        
        %make a mask for noiselevel
        
        subplot(2,ns,sl);
        imagesc(output.pars.xaxis,output.pars.yaxis,log10(abs(pic(:,:,sl))));
        set(gca,'Xticklabel','','Yticklabel','','clim',[2 4]);
        axis image;
        subplot(2,ns,ns+sl);
        imagesc(output.pars.xaxis,output.pars.yaxis,angle(pic(:,:,sl)));
        set(gca,'Xticklabel','','Yticklabel','');
        axis image;
        
        
    end
end


    