function output=deepi(timestamp,varargin)
% displacement-encoding-echo-planar-imaging
% out=deepi(timestamp,varargin)
% script to read the amplitude and phase of EPI images
% for the ligament/tendon project

options=varargin;

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
    time_flag(jj)=~isempty(strfind(fname,timestamp));
    img_flag(jj)=~isempty(strfind(fext,'img'));
    PH_flag(jj)=~isempty(strfind(fname,'_PH'));
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


output.image=amplitude.*phasefactor;
output.pars=outa.pars;
output.kspace=[];

    