%
% 
% A script to convert ginx format image files into a single analyze 
% format.
%
%
% A total rewrite of ginx2spm2 - after Tom Chenevert.
%
%
% Robert Welsh - Aug-17-00
% University of Michigan
% Department of Radiology
%  
%

% This scripts depends on : getginxvol
%                           getginx
%
%                           spm99/spm_hwrite.m
%                           anything that spm99/spm_hwrite uses.
%
%       (so in general this will work as long as you have spm99)
%


function ginx2spm3(saveFileName);

global ginxFOV
global ginxDIM

[ginxVol errmsg] = getginxvol;

ginxVOX = ginxFOV./ginxDIM;
ginxORG = ginxDIM / 2;
ginxSCL = 1;
ginxOFF = 0;
ginxDES = 'spm compatible - UMICH';

dataType = 4;

if errmsg == 1
  fprintf('Error in reading volume...');
else
  newVol = ginxVol;
  slices = size(ginxVol,3);
  for islice = 1:slices
    newVol(:,:,islice) = rot90(ginxVol(:,:,islice),2);
  end
  if nargin < 1
    saveFileName = input('Name to save file :');
  end
  if (isempty(findstr(saveFileName,'.img')))
    saveFileName = [saveFileName '.img'];  
  end
  fido = fopen(saveFileName,'w');
  count = fwrite(fido,newVol,'short');
  status = fclose(fido);
  spm_hwrite(saveFileName,ginxDIM,ginxVOX,ginxSCL,dataType,ginxOFF,ginxORG,ginxDES);
end

%
% All done.
%

