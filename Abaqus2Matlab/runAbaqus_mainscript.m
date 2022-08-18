clc;clear;

% File selection
% fileName = 'NoHoles_2_5MMDisp';
% fileName = 'NoHoles_5MMDisp';
% fileName = 'NoHoles_7MMDisp';
% fileName = '5MMHoles_2_5MMDisp';
fileName = '5MMHoles_5MMDisp';

addpath(genpath('Functions'))
addpath(genpath('InputFiles'))

[X,Y,Z,ElNode] = readInp(fileName);             % Step 1: Read Abaqus input
changeInp(fileName,length(X),length(ElNode))    % Step 2: Ammend input to include output requests
runInp(fileName)                                % Step 3: Run Abaqus through MatLab
addpath(genpath('Data'))
U = readDat(fileName,length(X));                % Step 4: Read Abaqus output data
save(['Data\' fileName '\' fileName '_test.mat'],'X','Y','Z','ElNode','U');
clc

% 08/18/2022
% To do: Make into function and feed into DefGrad_sim_v4