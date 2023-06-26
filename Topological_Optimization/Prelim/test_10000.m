clc;clear;
addpath(genpath('Functions'));
addpath(genpath('Data'));
load('test_settings_10000.mat');
output = sin_shear_opt(8,4,2,settings); 