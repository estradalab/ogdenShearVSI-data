%% Shape_optimization code
% Developed by Denislav P Nikolov (dnikolov@umich.edu)

clc;clear;
addpath(genpath('Functions'))
addpath(genpath('Data'))

% Name: Select a description for your code (for each time you run main)
desc_name = '10000el_optimMesh_x0_6_9_1.4';

% There are several options for function optimization to chose from
% func = 'test'; % Sample Rosenbrock function f(x) = 100*(x2-x1^2)^2 + (1-x1)^2
func = 'sin'; % Sine wave shape

% Optimization settings (parallel pool, display, and plot)
opts = optimset('fminsearch');
opts.UseParallel = true;
opts.Display = 'iter';
opts.PlotFcns = @optimplotfval;
opts.TolFun = 5e-2;

switch func
    case 'test'
        rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
        x0 = [0 0];
        LB = [2 -inf];
        UB = [inf 3];
        xsol = fminsearchbnd(rosen,x0,LB,UB,opts);
    case 'sin'
        % If you'd like to establish settings, open sin_shear_opt.m
        x0 = [6 9 1.4];
        LB = [6 0 0];
        UB = [20 10 8];
        A = [-1 0 2]; % Ensures that there is atleast 3mm thickness of material between sine wave valleys
        b = -3;
        tStart = tic;
        xsol = fminsearchcon(@ (x) sin_shear_opt_fun(x),x0,LB,UB,A,b,[],opts);
        xsol(1) = ceil(xsol(1));
        xsol(2) = round(xsol(2));
        xsol(3) = floor(xsol(3)*10^1)/10^1;
        tEnd = toc(tStart);
        optimization_folder = ['Data/' datestr(datetime(now,'ConvertFrom','datenum'),'yyyymmddTHHMMSS')];
        optimization_folder = [optimization_folder '_' desc_name];
        movefile('Data/Iter*',optimization_folder);
        save([optimization_folder '/optimization_settings.mat'],'x0','LB','UB','A','b','xsol','tEnd');
        saveas(gcf,[optimization_folder '/optimization_plot.png']); saveas(gcf,[optimization_folder '/optimization_plot.fig']);
end

% Calling sinusoidal geometry function
function f = sin_shear_opt_fun(x)
output = sin_shear_opt(ceil(x(1)),round(x(2)),floor(x(3)*10^1)/10^1);
% The single selected output is the inverse alpha goodness metric (trying to maximize)
f = 1/output.S(1);
end