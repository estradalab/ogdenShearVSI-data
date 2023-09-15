clc;clear;
warning off
addpath(genpath('Functions'));
addpath(genpath('Data'));

% test = 'test_100'; % Approximately 100 elements
% test = 'test_10000'; % Approximately 10000 elements
% test = 'treloar_test'; % Test if treloar data works for material parameters

% Neo-hookean parameters and lin/quad tets
% test = 'test_sweep_NH_nu_0.45_quad';
% test = 'test_sweep_NH_nu_0.25_quad';
% test = 'test_sweep_NH_nu_0.45_lin';
% test = 'test_sweep_NH_nu_0.25_lin';

% test = 'test_sweep_NH_uniaxial'; % Uniaxial extension study
% test = 'cost_function'; % Cost function space w/ fixed periods
% test = 'amp_sweep'; % Sweep between amplitudes
test = 'per_sweep'; % Sweep between periods
% test = 'element_test'; % Change number of elements for single trial

switch test
    case 'test_100'
        load('test_settings_100.mat');
        output = sin_shear_opt(8,4,0.6,settings);
    case 'test_10000'
        load('test_settings_10000.mat');
        output = sin_shear_opt(8,4,2,settings); 
    case 'treloar_test'
        load('test_settings_10000.mat');
        settings.save = 'treloar';
        settings.params = 'Arr-B-treloar';
        output = sin_shear_opt(8,2,2,settings);
        settings.params = 'NH-treloar';
        output = sin_shear_opt(8,2,2,settings);
        settings.params = 'M-R-treloar';
        output = sin_shear_opt(8,2,2,settings);
        movefile('Data/*treloar*',['Data/default_' num2str(settings.mesh_ref.num_of_el) 'el_Treloar']);
    case 'test_sweep_NH_uniaxial'
        load('test_settings_10000.mat');
        settings.mesh_ref.num_of_el = 10000;
        settings.params = 'neo-hooke-eco';
        settings.save = 'eco';
        settings.pres_disp = 8;
        settings.stretch = 'uniaxial';
        output = sin_shear_opt(8,0,0,settings);
        output = sin_shear_opt(8,2,0.25,settings);
        output = sin_shear_opt(8,2,0.5,settings);
        for i = 1
            for j = 1:3
                output = sin_shear_opt(8,2*j,i,settings);
            end
        end
        movefile('Data/Eco*',['Data/uniaxial_' num2str(settings.mesh_ref.num_of_el) 'el_Eco']);
    case {'test_sweep_NH_nu_0.45_quad','test_sweep_NH_nu_0.25_quad','test_sweep_NH_nu_0.45_lin','test_sweep_NH_nu_0.25_lin'}
        load('test_settings_10000.mat');
        settings.mesh_ref.num_of_el = 10000;
        switch test
            case {'test_sweep_NH_nu_0.45_quad','test_sweep_NH_nu_0.45_lin'}
                settings.params = 'neo-hooke-eco';
            case {'test_sweep_NH_nu_0.25_quad','test_sweep_NH_nu_0.25_lin'}
                settings.params = 'neo-hooke-eco-compr';
        end
        settings.save = 'eco';
        switch test
            case {'test_sweep_NH_nu_0.45_lin','test_sweep_NH_nu_0.25_lin'}
                settings.elementType = 'C3D4H';
        end
        output = sin_shear_opt(8,0,0,settings);
        output = sin_shear_opt(8,2,0.25,settings);
        output = sin_shear_opt(8,2,0.5,settings);
        output = sin_shear_opt(8,6,1,settings);
        for i = 1:2
            for j = 1:2
                output = sin_shear_opt(8,2*j,i,settings);
            end
        end
        movefile('Data/Eco*',['Data/' erase(test,'test_sweep_NH_') num2str(settings.mesh_ref.num_of_el) 'el_Eco']);
    case 'cost_function'
        load('test_settings_10000.mat');
        settings.save = 'delete';
        settings.params = 'neo-hooke-eco';

        %% Cost function space
        points = 100;
        d = linspace(6,10,floor(sqrt(points)));
        A = linspace(0,0.33,floor(sqrt(points)));
        S = zeros([length(d) length(A) 2]);
        tStart = tic;
        h = waitbar(0,'Progress: 0%');
        for i = 1:length(d)
            for j = 1:length(A)
                output = sin_shear_opt(d(i),3,A(j)*d(i),settings);
                S(i,j,1) = output.S(1);
                S(i,j,2) = output.S(2);
                T = seconds(toc(tStart));
                T.Format = 'hh:mm:ss';
                waitbar(2*((i-1)*length(d)+j)/numel(S),h,['Progress: ',num2str(floor(100*2*((i-1)*length(d)+j)/numel(S))),'% (Time Elapsed: ' char(T) ')'])
            end
        end
        tEnd = toc(tStart);
        close(h)
        save('Data/Cost_data/cost_function.mat','S','d','A','tEnd')
        
        %% Visualization
        contourf(d/settings.l,A,S(:,:,1)',10)
        colorbar
        xlabel('d/L [-], Normalized square cross section length','interpreter','latex')
        ylabel('A/d [-], Normalized sinusoidal amplitude','interpreter','latex')
        title('Goodness metric (Fixed sine wave @ 3 periods)','interpreter','latex')
        saveas(gcf,'Data/Cost_data/cost_function.png'); saveas(gcf,'Data/Cost_data/cost_function.pdf')

    case 'per_sweep'
        load('test_settings_10000.mat');
        settings.save = 'delete';
        settings.params = 'neo-hooke-eco';

        sin_sweep = 0:10;
        tStart = tic;
        h = waitbar(0,'Progress: 0%');
        for i = 1:length(sin_sweep)
            output = sin_shear_opt(6,sin_sweep(i),2,settings);
            S(i,1) = output.S(1);
            S(i,2) = output.S(2);
            T = seconds(toc(tStart));
            T.Format = 'hh:mm:ss';
            waitbar(i/length(sin_sweep),h,['Progress: ',num2str(floor(100*i/length(sin_sweep))),'% (Time Elapsed: ' char(T) ')'])
        end
        tEnd = toc(tStart);
        close(h)
        save('Data/Sin_data/per_sweep.mat','S','sin_sweep')
        plot(sin_sweep,S(:,1))
        xlabel('N [-], Sinusoidal periods','interpreter','latex')
        ylabel('S [-], Normalized goodness metric','interpreter','latex')
        title('Goodness metric (Varied sine wave @ A = 2 mm, d = 6 mm)','interpreter','latex')
        saveas(gcf,'Data/Sin_data/per_sweep.png'); saveas(gcf,'Data/Sin_data/per_sweep.pdf')
    case 'amp_sweep'
        load('test_settings_10000.mat');
        settings.save = 'delete';
        settings.params = 'neo-hooke-eco';
        

        sin_sweep = linspace(0,0.33,100);
        tStart = tic;
        h = waitbar(0,'Progress: 0%');
        for i = 1:length(sin_sweep)
            output = sin_shear_opt(6,3,sin_sweep(i)*6,settings);
            S(i,1) = output.S(1);
            S(i,2) = output.S(2);
            T = seconds(toc(tStart));
            T.Format = 'hh:mm:ss';
            waitbar(i/length(sin_sweep),h,['Progress: ',num2str(floor(100*i/length(sin_sweep))),'% (Time Elapsed: ' char(T) ')'])
        end
        tEnd = toc(tStart);
        close(h)
        save('Data/Cost_data/sin_sweep.mat','S','sin_sweep')
        plot(sin_sweep,S(:,1))
        xlabel('A/d [-], Normalized sinusoidal amplitude','interpreter','latex')
        ylabel('S [-], Normalized goodness metric','interpreter','latex')
        title('Goodness metric (Fixed sine wave @ 3 periods, d = 6 mm)','interpreter','latex')
        saveas(gcf,'Data/Cost_data/amplitude_sweep.png'); saveas(gcf,'Data/Cost_data/amplitude_sweep.pdf')
    case 'element_test'
        load('test_settings_10000.mat');
        settings.save = 'delete';
        settings.params = 'neo-hooke-eco';
        settings.save = 'eco';
        
        el_sweep = [100 200 500 1000 1500 2000 5000 10000 12000 15000 20000];
        for i = 1:length(el_sweep)
            settings.mesh_ref.num_of_el = el_sweep(i);
            output = sin_shear_opt(6,3,2,settings);
            S(i,1) = output.S(1);
            S(i,2) = output.S(2);
        end
        save('Data/Cost_data/element_sweep.mat','S','el_sweep')
        plot(el_sweep,S(:,1))
end