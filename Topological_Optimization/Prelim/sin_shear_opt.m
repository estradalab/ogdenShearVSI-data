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
% X_ref - [X,Y,Z] coordinates of all nodes (X_ref.node) and elements (X_ref.el)
%       reference position array is in the form X_ref.[node/el](node,dir_[1/2/3])
% U - [U1,U2,U3] displacements in all three cartesian direction of all nodes
%       displacement array is in the form U(node,dir_[1/2/3])
% RF1 - reaction force along boundary controlled surfaces
%       RF1.bc1 is displaced surface; RF1.bc2 is encastered surface
%       All other elements that are not of the boundary surfaces are NaN values
% traction - traction required to pull specimen in shear (units - Newtons)
% F_t - deformation gradient for a single time step for all elements
%       deformation gradient is in the form F_t{time}{i,j}(element)
% k/k3 - stretch state decomposition for all elements (custom decomp/criscione/decomp)
% lam/k2 - maximum stretch amplitude decomposition for all elements (custom decomp/criscione/decomp)
%       decomposition parameters are in the form k/lam{time}(element)
% Salpha_all - alpha sensitivity map for pre-selected alpha parameter
% Smu_all - mu sensitivity map for pre-selected alpha parameter
%       sensitivity maps are in the form S_all(lambda,alpha,k)
% H - 2D histogram of simulation (custom decomp only, for criscione [in developement])
%       histograms are in the form H(lambda,k)
% S(1) - alpha goodness metric; S(2) - mu goodness metric

% Computer specs for approximate times:
%   Processor: Intel(R) Xeon(R) CPU E5-1620 v2 @ 3.70GHz
%   Ram: 16.0 GB
%   OS: Windows 10 Enterprise
%   Abaqus: 2021 standard
%   MATLAB: R2022a

% For test case in the approximate time calculations, use 
% Run test_100.m

% For fast test case, use
% Run test_10000.m

% If you'd like to rerun test cases, files will be overwritten, so do save
% locally before running multiple times.

% Set parameters
if isempty(varargin)
    % Current settings are meant for the Nelder-Mead optimizer
    settings.params = 'neo-hooke-eco'; % Utilizing treloar data for ogden model parameters
    settings.mesh_ref.num_of_el = 10000; % Sets approximate number of elements of the mesh
    settings.pres_disp = 3; % Prescribed displacement of 3 mm
    settings.mesh = 'tet'; % Quadratic tetrahedral elements (other element types are in progress)
    settings.l = 40; % Sample length is 40 mm
    settings.abaqus_ver = '2021'; % Abaqus version
    settings.elementType = 'C3D10H'; % Element type: C3D10H is quadratic hybrid tet elements
    settings.alpha = 2; % Alpha parameter of Ogden model used to calculate sensitivity mapping
    settings.maxLam = 2; % Maximum lambda used for the 2D histogram and sensitivity plots (2 decimal point-limit)
    settings.bin_res = 0.01; % Bin resolution for 2D histogram and sensitivity plot
    settings.parallel = false; % Uses parfor loop in decomposition loop. Required as false for optimization
    settings.sigma_calc = false; % Optional acquisition of sigma tensor for entire sample (outdated)
    settings.save = 'optim'; % 'none' - doesn't save data; 'test' - saves data for test; 'optim' - saves data for optimization runs
    settings.mesh_ref.exact = false; % Iterates size of mesh element until number of elements is exact (Works for 10000)
    settings.stretch = 'shear'; % 'shear' is to pull specimen in shear on the sinusoidal profile, and 'uniaxial' is to pull it in uniaxial extension
else
    settings = varargin{1};
end

fileName = ['sq-' num2str(d) 'mm_sin-per-' num2str(N) '_sin-amp-' num2str(A) 'mm_' settings.mesh];

% Approximate processing time for test settings mesh generation: 23 seconds
[X,Y,Z,ElNode] = genMesh(fileName,settings.mesh,settings.params,settings.pres_disp,...
    settings.mesh_ref,settings.abaqus_ver,settings.elementType,d,N,A,settings.l,settings.stretch);
output.X_ref.node = [X,Y,Z]; % Nodal positions
% Element centroidal positions
for idx = 1:length(ElNode)
    for i = 1:8
        n(i) = ElNode(idx,i);
    end
    output.X_ref.el(idx,:) = [mean(X(n(:))),mean(Y(n(:))),mean(Z(n(:)))];
end

% Approximate processing time for test settings simulation: 99 seconds
% Output deformation gradient
[output.U,F,output.RF1] = runAbaqus(fileName,X);
movefile([fileName '.inp'],['Data/' fileName]);

% Calculation for traction force
output.traction = nansum(output.RF1.bc1);

% Approximate processing time for test settings decomposition calculation: 1 seconds
% Carry out the decomposition of F to k/lambdabda and k3/k2
output.F_t{1} = F;
[output.k,output.lam,~] = param_decoup_main(output.F_t,settings.parallel,'ftolamandk');
[output.k3,output.k2,~] = param_decoup_main(output.F_t,settings.parallel,'ftok2andk3');

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

switch settings.save
    case 'none'
    case 'delete'
        rmdir 'Data/sq*' s
    case 'test'
        writematrix(output.X_ref.node,['Data/' fileName '/x_ref_' fileName '.csv']); 
        writematrix(output.U,['Data/' fileName '/node_disp_' fileName '.csv']); 
        save(['Data/' fileName '/output.mat'],'output');
        movefile(['Data/' fileName],['Data/Test_' fileName]);
    case 'optim'
        writematrix(output.X_ref.node,['Data/' fileName '/x_ref_' fileName '.csv']); 
        writematrix(output.U,['Data/' fileName '/node_disp_' fileName '.csv']); 
        save(['Data/' fileName '/output.mat'],'output');
        movefile(['Data/' fileName],['Data/Iter' num2str(length(dir('Data/Iter*'))+1) '_' fileName]);
    case 'eco'
        writematrix(output.X_ref.node,['Data/' fileName '/x_ref_' fileName '.csv']); 
        writematrix(output.U,['Data/' fileName '/node_disp_' fileName '.csv']); 
        save(['Data/' fileName '/output.mat'],'output');
        movefile(['Data/' fileName],['Data/Eco_' fileName]);
end

end