% test = 'test_100'; % Approximately 100 elements
% test = 'test_10000'; % Approximately 10000 elements
test = 'test_sweep_NH'; % Neo-hookean parameters
% test = 'cost_function'; % Cost function space w/ fixed periods
% test = 'sin_sweep'; % Sweep between periods
% test = 'element_test'; % Change number of elements for single trial

switch test
    case 'test_100'
        clc;clear;
        warning off
        addpath(genpath('Functions'));
        addpath(genpath('Data'));
        load('test_settings_100.mat');
        output = sin_shear_opt(8,4,0.6,settings);
    case 'test_10000'
        clc;clear;
        warning off
        addpath(genpath('Functions'));
        addpath(genpath('Data'));
        load('test_settings_10000.mat');
        output = sin_shear_opt(8,4,2,settings); 
    case 'test_sweep_NH'
        clc;clear;
        warning off
        addpath(genpath('Functions'));
        addpath(genpath('Data'));
        load('test_settings_10000.mat');
        settings.params = 'neo-hooke-eco';
        settings.save = 'eco';
        output = sin_shear_opt(8,0,0,settings);
        for i = 1:3
            for j = 1:2
                output = sin_shear_opt(8,2*j,i,settings);
            end
        end
    case 'cost_function'
        clc;clear;
        warning off
        addpath(genpath('Functions'));
        addpath(genpath('Data'));
        load('test_settings_10000.mat');
        settings.save = 'delete';

        %% Cost function space
        points = 1000;
        d = linspace(6,10,floor(sqrt(points)));
        A = linspace(0,1,floor(sqrt(points)));
        S = zeros([length(d) length(A) 2]);
        tStart = tic;
        h = waitbar(0,'Progress: 0%');
        for i = 1:length(d)
            for j = 1:length(A)
                output = sin_shear_opt(d(i),3,A(j),settings);
                S(i,j,1) = output.S(1);
                S(i,j,2) = output.S(2);
                T = seconds(toc(tStart));
                T.Format = 'hh:mm:ss';
                waitbar(2*((i-1)*length(d)+j)/numel(S),h,['Progress: ',num2str(floor(100*2*((i-1)*length(d)+j)/numel(S))),'% (Time Elapsed: ' char(T) ')'])
            end
        end
        tEnd = toc(tStart);
        close(h)
        save('Data/test_data/cost_function.mat','S','d','A','tEnd')
        
        %% Visualization
        contourf(A,d,S(:,:,1))
        colorbar
        ylabel('d [mm], Square cross section length','interpreter','latex')
        xlabel('A [mm], Sinusoidal amplitude','interpreter','latex')
        title('Goodness metric (Fixed sine wave @ 3 periods)','interpreter','latex')
        saveas(gcf,'Data/test_data/cost_function.png'); saveas(gcf,'Data/test_data/cost_function.pdf')
    case 'sin_sweep'
        clc;clear;
        warning off
        addpath(genpath('Functions'));
        addpath(genpath('Data'));
        load('test_settings_10000.mat');
        settings.save = 'delete';

        sin_sweep = 0:10;
        clear S
        for i = sin_sweep
            output = sin_shear_opt(6,i,2,settings);
            S(i+1,1) = output.S(1);
            S(i+1,2) = output.S(2);
        end
        save('Data/test_data/sin_sweep.mat','S','sin_sweep')
        plot(sin_sweep,S(:,1))
    case 'element_test'
        clc;clear;
        warning off
        addpath(genpath('Functions'));
        addpath(genpath('Data'));
        load('test_settings_10000.mat');
        settings.save = 'delete';
        
        el_sweep = [100 200 500 1000 1500 2000 5000 10000 12000 15000 20000];
        for i = 1:length(el_sweep)
            settings.mesh_ref.num_of_el = el_sweep(i);
            output = sin_shear_opt(6,3,2,settings);
            S(i,1) = output.S(1);
            S(i,2) = output.S(2);
        end
        save('Data/test_data/sin_sweep.mat','S','el_sweep')
        plot(el_sweep,S(:,1))
end