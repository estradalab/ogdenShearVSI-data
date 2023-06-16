function [U,F] = runAbaqus(fileName,varargin)

% Step 1: Ammend input to include output requests
X = varargin{1}; Y = varargin{2}; Z = varargin{3}; ElNode = varargin{4};
changeInp(fileName,length(X),length(ElNode))

% Step 2: Run Abaqus through MatLab
% Approximate time to run 40,000 element mesh simulation: 41 seconds
runInp(fileName)

% Step 3: Read Abaqus output data
addpath(genpath('Data'))
[U,F] = readDat(fileName,length(X));
clc