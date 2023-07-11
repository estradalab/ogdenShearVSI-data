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
% S12 - shear stress along boundary controlled surfaces
%       S12.bc1 is displaced surface; S12.bc2 is encastered surface
%       All other elements that are not of the boundary surfaces are NaN values
% sig - stress tensor for all elements
%       stress tensor is in the form sig{i,j}(element)
% traction - traction required to pull specimen in shear (units - Newtons)
% F_t - deformation gradient for a single time step for all elements
%       deformation gradient is in the form F_t{time}{i,j}(element)
% k - stretch state decomposition for all elements
% lam - maximum stretch amplitude decomposition for all elements
%       decomposition parameters are in the form k/lam{time}(element)
% Salpha_all - alpha sensitivity map for pre-selected alpha parameter
% Smu_all - mu sensitivity map for pre-selected alpha parameter
%       sensitivity maps are in the form S_all(lambda,alpha,k)
% H - 2D histogram of simulation
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
    settings.params = 'ogden-treloar'; % Utilizing treloar data for ogden model parameters
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
    settings.sigma_calc = false; % Optional acquisition of sigma tensor for entire sample
    settings.save = 'optim'; % 'none' - doesn't save data; 'test' - saves data for test; 'optim' - saves data for optimization runs
    settings.mesh_ref.exact = true; % Iterates size of mesh element until number of elements is exact (Works for 10000)
else
    settings = varargin{1};
end

fileName = ['sq-' num2str(d) 'mm_sin-per-' num2str(N) '_sin-amp-' num2str(A) 'mm_' settings.mesh];

% Approximate processing time for test settings mesh generation: 23 seconds
[X,Y,Z,ElNode] = genMesh(fileName,settings.mesh,settings.params,settings.pres_disp,...
    settings.mesh_ref,settings.abaqus_ver,settings.elementType,d,N,A,settings.l);
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
[output.U,F,output.S12,output.sig] = runAbaqus(fileName,X,settings.sigma_calc);
movefile([fileName '.inp'],['Data/' fileName]);

% Calculation for traction force
% First, calculate arc length of sine wave
x = linspace(0,settings.l,100); dx = diff(x);
y = A*sin(2*pi*N*x/settings.l); dy = diff(y);
arc_l = sum(sqrt(dx.^2+dy.^2));
% Then, backward calculate V = tau*A to determine shear force, and thus traction
area = arc_l*d;
output.traction = area*nanmean([nanmean(output.S12.bc1) nanmean(output.S12.bc2)]);

% Approximate processing time for test settings decomposition calculation: 1 seconds
% Carry out the decomposition of F to k/lambdabda
output.F_t{1} = F;
[output.k,output.lam,~] = param_decoup_main(output.F_t,settings.parallel);

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
end

end