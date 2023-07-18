function [U,F,RF1] = runAbaqus(fileName,varargin)

% Step 1: Ammend input to include output requests
X = varargin{1};
changeInp(fileName)

% Step 2: Run Abaqus through MatLab
% Approximate time to run 40,000 element mesh simulation: 41 seconds
runInp(fileName)

% Step 3: Read Abaqus output data
addpath(genpath('Data'))
[U,F,RF1] = readDat(fileName,length(X));