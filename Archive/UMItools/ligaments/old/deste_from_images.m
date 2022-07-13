function output=deste_from_images(timestamp,varargin)
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

%get the amplitude;
outa=varianms(img_directory,options{:});
amplitude=outa.image;

%get the phase;
outp=varianms(phi_directory,options{:});
phasefactor=exp(1i*outp.image);

%get the paramters out of the fid directory
pars=getparameters([fid_directory '/procpar']);

%construct the complex image
output.image=amplitude.*phasefactor;

%extract the required parameters for plotting, construct the axes
%all lengths in cm
[numro numpe numsl]=size(amplitude);
output.pars.xaxis=[0:numpe-1]/numpe*pars.lpe;
output.pars.yaxis=[0:numro-1]/numro*pars.lro;
output.pars.zaxis=pars.pss;
output.pars.lambda=pars.lambda;

%claculate the diffusive attenuation factor b=k^2 Delta
%standard units are s/mm^2
output.pars.diffusive_b=(2*pi/(pars.lambda*10))^2*pars.TS;


    