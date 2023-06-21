clc;clear;
addpath(genpath('Functions'));
load('test_settings_100.mat');
output = sin_shear_opt(8,4,0.6,settings); 
fileName = 'sq-8mm_sin-per-4_sin-amp-0.6mm_tet'; 
writematrix(output.X_ref.node,['Data/' fileName '/x_ref_' fileName '.csv']); 
writematrix(output.U,['Data/' fileName '/node_disp_' fileName '.csv']); 
save(['Data/' fileName '/output.mat'],'output');
