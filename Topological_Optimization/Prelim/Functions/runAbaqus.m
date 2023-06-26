function [U,F,S12,sig] = runAbaqus(fileName,varargin)

% Step 1: Ammend input to include output requests
X = varargin{1}; sigma_calc = varargin{2};
changeInp(fileName,sigma_calc)

% Step 2: Run Abaqus through MatLab
% Approximate time to run 40,000 element mesh simulation: 41 seconds
runInp(fileName)

% Step 3: Read Abaqus output data
addpath(genpath('Data'))
[U,F,S12,sig] = readDat(fileName,length(X),sigma_calc);