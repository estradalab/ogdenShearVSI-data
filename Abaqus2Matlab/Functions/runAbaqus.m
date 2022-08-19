function [X,Y,Z,ElNode,U] = runAbaqus(fileName)

[X,Y,Z,ElNode] = readInp(fileName);                 % Step 1: Read Abaqus input
if ~exist(['Data\' fileName '\' fileName '_test.dat'])
    changeInp(fileName,length(X),length(ElNode))    % Step 2: Ammend input to include output requests
    runInp(fileName)                                % Step 3: Run Abaqus through MatLab
end
addpath(genpath('Data'))
U = readDat(fileName,length(X));                    % Step 4: Read Abaqus output data
clc