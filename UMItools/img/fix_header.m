%
% R.C.Welsh / Univ of Michigan
% Department of Radiology
% May 11 - 2000
%
%
% A little matlab script to take in a header and rewrite it out
% using the correct type. 
%
% This is to fix a problem with analyze headers that are written by 
% us. We put out data type 1024, but then AIR thinks this is 4096
% bits/pixel. So I will use spm_hread and spm_hwrite to fix the
% problem.
%
% Basically - 1024 -> int16 -> 4 is what needs to be done.
%
% v.01 - Do just the header rebuild.
%   
% v.02 - Nope, header just screws up little/big endian.
%        Trick now is just to open and read the file and 
%        write the file. SPM/MatLab will handle the conversion 
%        Between little and big endian
%
%        This little script was tested and AIR seemd to like it.
%
%               a = spm_vol('s1_oct26_99_e5501s8r1_007') 
%               [y xyz]= spm_read_vols(a));
%               spm_write_vol(a,y)

%
% use spm_get to point to the files.
% 

FileList = spm_get([0 9999],'img','Pick images to fix','./');

numFiles = size(FileList);
numFiles = numFiles(1);

fprintf('Fixing files ');

for iFile = 1:numFiles
  FileName              = FileList(iFile,:);
  aHeader               = spm_vol(FileName);
  [aYMatrix aXYZMatrix] = spm_read_vols(aHeader);
  spm_write_vol(aHeader,aYMatrix);
  fprintf('.');
end

fprintf(' done\n');

%
% All done.
%
