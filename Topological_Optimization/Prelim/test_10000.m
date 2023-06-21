clc;clear;
addpath(genpath('Functions'));
load('test_settings_10000.mat');
output = sin_shear_opt(8,4,2,settings); 
fileName = 'sq-8mm_sin-per-4_sin-amp-2mm_tet';
writematrix(output.X_ref.node,['Data/' fileName '/x_ref_' fileName '.csv']); 
writematrix(output.U,['Data/' fileName '/node_disp_' fileName '.csv']); 
save(['Data/' fileName '/output.mat'],'output');