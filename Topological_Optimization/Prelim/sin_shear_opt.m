function output = sin_shear_opt(d,N,A,varargin)
% Function that calculates goodness metric (S) at a select Ogden parameter.
% Note that all length parameters are in mm.

% Inputs:
% d - square cross-section dimension (mm)
% N - sinusoidal periods in sample
% A - sinusoidal amplitude of sample (mm)

% Optional parameter (you can input your own parameters below)
% varargin - test settings (provided as test_settings.mat as provided)

% Output structure:
% X_ref_node - [X,Y,Z] coordinates of all nodes
%       reference position array is in the form X_ref_node(node,dir_[1/2/3])
% U - [U1,U2,U3] displacements in all three cartesian direction of all nodes
%       displacement array is in the form U(node,dir_[1/2/3])
% F_t - deformation gradient for a single time step for all elements
%       deformation gradient is in the form F_t{time}{i,j}(element)
% k - stretch state decomposition for all elements
% lam - maximum stretch amplitude decomposition for all elements
%       decomposition parameters are in the form k/lam{time}(element)
% Salpha_all - alpha sensitivity map for pre-selected alpha parameter
% Smu_all - mu sensitivity map for pre-selected alpha parameter
%       sensitivity maps are in the form S_all(lambda,alpha,k)

% Computer specs for approximate times:
%   Processor: Intel(R) Xeon(R) CPU E5-1620 v2 @ 3.70GHz
%   Ram: 16.0 GB
%   OS: Windows 10 Enterprise
%   Abaqus: 2021 standard
%   MATLAB: R2022a

% For test case in the approximate time calculations, use 
% clc;clear;addpath(genpath('Functions'));load('test_settings_10000.mat');output = sin_shear_opt(8,4,2,settings); fileName = 'sq-8mm_sin-per-4_sin-amp-2mm_tet'; save(['Data/' fileName '/output.mat'],'output');

% For fast test case, use
% clc;clear;addpath(genpath('Functions'));load('test_settings_100.mat');output = sin_shear_opt(8,4,0.6,settings); fileName = 'sq-8mm_sin-per-4_sin-amp-0.6mm_tet'; save(['Data/' fileName '/output.mat'],'output');

% If you'd like to rerun test cases, make sure to move or delete the files
% rmdir(['Data/' fileName], 's')

addpath(genpath('Functions'))

% Set parameters
if isempty(varargin)
    settings.params = 'ogden-treloar'; % Utilizing treloar data for ogden model parameters
    settings.mesh_ref.num_of_el = 10000; % Sets approximate number of elements of the mesh
    settings.pres_disp = 3; % Prescribed displacement of 3 mm (minimum of a = 8 mm is required)
    settings.mesh = 'tet'; % Quadratic tetrahedral elements (other element types are in progress)
    settings.l = 40; % Sample length is 40 mm
    settings.abaqus_ver = '2021'; % Abaqus version
    settings.elementType = 'C3D10H'; % Element type: C3D10H is quadratic hybrid tet elements
    settings.alpha = 2; % Alpha parameter of Ogden model used to calculate sensitivity mapping
    settings.maxLam = 2; % Maximum lambda used for the 2D histogram and sensitivity plots (2 decimal point-limit)
    settings.bin_res = 0.01; % Bin resolution for 2D histogram and sensitivity plot
else
    settings = varargin{1};
end

fileName = ['sq-' num2str(d) 'mm_sin-per-' num2str(N) '_sin-amp-' num2str(A) 'mm_' settings.mesh];

% Approximate processing time for test settings mesh generation: 10 seconds
[X,Y,Z,ElNode] = genMesh(fileName,settings.mesh,settings.params,settings.pres_disp,...
    settings.mesh_ref,settings.abaqus_ver,settings.elementType,d,N,A,settings.l);
output.X_ref_node = [X,Y,Z]; % Nodal positions
% Approximate processing time for test settings simulation: 150 seconds
% Output deformation gradient
[output.U,F] = runAbaqus(fileName,X,Y,Z,ElNode);
movefile([fileName '.inp'],['Data/' fileName]);

% Approximate processing time for test settings decomposition calculation: 4 seconds
% Carry out the decomposition of F to k/lambdabda
output.F_t{1} = F;
[output.k,output.lam,~] = param_decoup_main(output.F_t);

% Store the lambda as a column vector
lambda = 1:settings.bin_res:settings.maxLam;

% Sensitivity functions with respect to alpha and mu
mu = 1; % Normalizing mu
Salpha_gen = @(k) mu*bsxfun(@times,log(lambda'),((-k^2+1.5*k+0.5)^2*bsxfun(@power,lambda',(-k^2+1.5*k+0.5)*settings.alpha-1)...
    +(-k+0.5)^2*bsxfun(@power,lambda',(-k+0.5)*settings.alpha-1)...
    +(k^2-0.5*k-1)^2*bsxfun(@power,lambda',(k^2-0.5*k-1)*settings.alpha-1)));
Smu_gen = @(k) (-k^2+1.5*k+0.5)*bsxfun(@power,lambda',(-k^2+1.5*k+0.5)*settings.alpha-1)...
    +(-k+0.5)*bsxfun(@power,lambda',(-k+0.5)*settings.alpha-1)...
    +(k^2-0.5*k-1)*bsxfun(@power,lambda',(k^2-0.5*k-1)*settings.alpha-1);

% Create 2D array with the senstivity functions given k and lambda range
k_=0:settings.bin_res:1;
for j = 1:length(k_)
    output.Salpha_all(:,:,j) = Salpha_gen(k_(j));
    output.Smu_all(:,:,j) = Smu_gen(k_(j));
end

% Make 2D histogram data
[~,~,h] = ndhist(output.lam{1}(:),output.k{1}(:),'axis',[1,0;settings.maxLam,1],'fixed_bin_res',settings.bin_res);
output.H = h';

% Cross-convolution to determine goodness metric, S(1) = S_alpha; S(2) = S_mu;
mat_temp = squeeze(output.Salpha_all(:,1,:));
mat_temp2 = squeeze(output.Smu_all(:,1,:));
nrmlz = sum(output.H(:));
output.S(1) = sum(output.H(:).*mat_temp(:))/nrmlz;
output.S(2) = sum(output.H(:).*mat_temp2(:))/nrmlz;

end