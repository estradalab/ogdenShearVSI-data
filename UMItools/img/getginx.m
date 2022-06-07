%
% A little function to read an ginx file.
% 
% Some code due to Tom Chenevert.

% Robert Welsh - Aug-17-00
% University of Michigan
% Department of Radiology
%  
%

%
% format:
%
function [yImg,numimages,dim1,dim2,errmsg] = getginx(imageFileName)

global ginxFOV
global ginxDIM
global ginxVOX

fid = fopen(imageFileName,'r','b');

errmsg = 0;

if ( fid > 0 ) 
  numim_off = 1368; % 1368 Offset to num of images in series
  status = fseek(fid,numim_off,'bof'); % Byte offset from current file position.
  numimages = fread(fid,1,'long'); 
  % Total number of images in scan-push for this series
  %numimages = 54; % Temphardcode 54 to avoid DW xyz images
  
  exno_off = 2180; % 2180+12 Offset to exam number
  status = fseek(fid,exno_off+12,'bof'); 
  % Byte offset from current file position.
  exno = fread(fid,1,'short'); %%
  serno = fread(fid,1,'short');%%
  imno = fread(fid,1,'short');%%
  slthk_off = exno_off+32; % 2180+12 Offset to exam number
  status = fseek(fid,slthk_off,'bof'); % Byte offset from current file position.
  slthk = fread(fid,1,'float');%%
  matx = fread(fid,1,'short');%%; % Is recontructed/image matrix size
  maty = fread(fid,1,'short');%%; % Is recontructed/image matrix size
  dfovx = fread(fid,1,'float');%%;
  dfovy = fread(fid,1,'float');%%;
  dimx = fread(fid,1,'float');%%; % Is acquired matrix size
  dimy = fread(fid,1,'float');%%; % Is acquired matrix size
  % matrix=input('Matrix size of images 128 or 256 ?');
  dim1 = matx;
  dim2 = maty;
  nimg = numimages; %dim3;
  status = fseek(fid,-2*dim1*dim2,1);
  yImg = fread(fid,[dim1,dim2],'short');
  ginxFOV = [dfovx dfovy numimages*slthk];
  ginxDIM = [dim1 dim2 numimages];
  fclose(fid);
else
  errmsg = 1;
end

%
% Done
%
