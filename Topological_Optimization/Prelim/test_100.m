clc;clear;
warning off
addpath(genpath('Functions'));
addpath(genpath('Data'));
load('test_settings_100.mat');
output = sin_shear_opt(8,4,0.6,settings); 
