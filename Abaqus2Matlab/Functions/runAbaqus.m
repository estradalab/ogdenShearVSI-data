function [X,Y,Z,ElNode,U] = runAbaqus(curdir,fileName,varargin)

if ~isempty(varargin)
    X = varargin{1}; Y = varargin{2}; Z = varargin{3}; ElNode = varargin{4};
else
    if contains(curdir,'sweep','IgnoreCase',true)                                                   % Step 1: Read Abaqus input
        [X,Y,Z,ElNode] = readInp([curdir '\' fileName]);
    else
        [X,Y,Z,ElNode] = readInp(fileName);
    end
end

if contains(curdir,'sweep','IgnoreCase',true)
    if ~exist(['Data\' curdir '\' fileName '\' fileName '_test.dat'])
        changeInp([curdir '\' fileName],fileName,length(X),length(ElNode))  % Step 2: Ammend input to include output requests
        runInp([curdir '\' fileName],fileName)                                                      % Step 3: Run Abaqus through MatLab
    end
else
    if ~exist(['Data\' fileName '\' fileName '_test.dat'])
        changeInp(fileName,fileName,length(X),length(ElNode))
        runInp(fileName,fileName)
    end
end

addpath(genpath('Data'))
U = readDat(fileName,length(X));                                                                        % Step 4: Read Abaqus output data
clc