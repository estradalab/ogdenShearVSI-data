function [X,Y,Z,ElNode,U] = runAbaqus(fileName,varargin)

if ~isempty(varargin)
    X = varargin{1}; Y = varargin{2}; Z = varargin{3}; ElNode = varargin{4};
else
    [X,Y,Z,ElNode] = readInp(fileName);                 % Step 1: Read Abaqus input
end

if ~exist(['Data\' fileName '\' fileName '_test.dat'])
    changeInp(fileName,length(X),length(ElNode))    % Step 2: Ammend input to include output requests
    runInp(fileName)                                % Step 3: Run Abaqus through MatLab
end
addpath(genpath('Data'))
U = readDat(fileName,length(X));                    % Step 4: Read Abaqus output data
clc