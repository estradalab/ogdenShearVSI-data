clc; clear;

addpath(genpath('Functions'))
addpath(genpath('InputFiles'))

% curdir = '22-0201-5MM_Holes_SS';
% curdir = '22-0301-NoHoles_SS';
% curdir = '22-0325-Uniaxial';
% curdir = '22-1212-Wavy_SS';
% curdir = '22-1215-Wavy_Sweep';
% curdir = '23-0119-Wavy_Prd_Sweep';
% curdir = '23-0202-Wavy-Thick_Amp_Sweep';
% curdir = '23-0202-Wavy-Thick_Prd_Sweep';
% curdir = '23-0217-Wavy-Thick_Width_Sweep';
curdir = '23-0405-SingleCompression';

params = 'ogden-treloar';

mesh_ref.defsize = 0.5; % 50% of default size elements in mesh
mesh_ref.maxelsize = 0.3; % Keeps a constant element density throughout the sweep  (overrides mesh_ref.defsize; comment out if defsize is preferred)
% Ran mesh_ref.maxelsize for 0217 and later

saveset = true;
generate_mesh = true; % Generate mesh (for false, you'll be using hexahedral elements pre-made in Abaqus)

switch curdir
    case '22-0201-5MM_Holes_SS'
        A{1} = '5MMHoles_2_5MMDisp';
        A{2} = '5MMHoles_5MMDisp';
        pres_disp{1} = 2.5; pres_disp{2} = 5;
    case '22-0301-NoHoles_SS'
        A{1} = 'NoHoles_2_5MMDisp';
        A{2} = 'NoHoles_5MMDisp';
        A{3} = 'NoHoles_7MMDisp';
        pres_disp{1} = 2.5; pres_disp{2} = 5; pres_disp{3} = 7;
    case '22-0325-Uniaxial'
        A{1} = 'Uniaxial_7MMDisp';
        pres_disp{1} = 7;
    case '22-1212-Wavy_SS'
        A{1} = 'ShearRect_6_25MMDisp';
        A{2} = 'ShearWavy_6_25MMDisp';
        pres_disp{1} = 6.25; pres_disp{2} = 6.25;
    case '22-1215-Wavy_Sweep'
         sin_sweep = linspace(0,2,11);
        for iii = 1:length(sin_sweep)
            A{iii} = ['ShearWavy_6.25MMDisp_Amp_' num2str(sin_sweep(iii))];
            pres_disp{iii} = 6.25;
        end
    case '23-0119-Wavy_Prd_Sweep'
        sin_sweep = 1:10;
        for iii = 1:length(sin_sweep)
            A{iii} = ['ShearWavyPrd1_6.25MMDisp_Amp_0.6MM_Prd_' num2str(sin_sweep(iii))];
            pres_disp{iii} = 6.25;
        end
    case '23-0202-Wavy-Thick_Amp_Sweep'
        sin_sweep = linspace(0,2,11);
        for iii = 1:length(sin_sweep)
            A{iii} = ['ShearWavyAmp2_6.25MMDisp_Amp_' num2str(sin_sweep(iii))];
            pres_disp{iii} = 6.25;
        end
    case '23-0202-Wavy-Thick_Prd_Sweep'
        sin_sweep = 1:10;
        for iii = 1:length(sin_sweep)
            A{iii} = ['ShearWavyPrd2_6.25MMDisp_Amp_0.6MM_Prd_' num2str(sin_sweep(iii))];
            pres_disp{iii} = 6.25;
        end
    case '23-0217-Wavy-Thick_Width_Sweep'
        width_sweep = 5:5:60;
        for iii = 1:length(width_sweep)
            A{iii} = ['ShearWavyWidth1_6.25MMDisp_Amp_0.6MM_Prd_4_Width_' num2str(width_sweep(iii))];
            pres_disp{iii} = 6.25;
        end
    case '23-0405-SingleCompression'
        A{1} = ['SingleCompression_6X6X1_0.15MMDisp'];
        pres_disp{1} = -0.15;
    case '22-0516-Ogden_Methodical'
        % Work in progress
    case '22-0526-Ogden_Methodical_Size'
        % Work in progress
    case '22-0615-Ogden_Iosipescu'
        % Work in progress
    case '22-0623-Ogden_Iosipescu_Magnitude'
        % Work in progress
    case '22-0623-Ogden_Iosipescu_HalfConstr'
        % Work in progress
end

if generate_mesh
    mesh = 'tet';
    for i = 1:length(A)
        if ~(isfile(['InputFiles\' A{i} '_tet.inp']) || isfile(['InputFiles\' curdir '\' A{i} '_tet.inp']))
            [x_temp,y_temp,z_temp,TRI_temp] = genMesh(curdir,[A{i} '_tet'],mesh,params,'Write',pres_disp{i},mesh_ref);
        else
            switch A{i}
                case {'ShearWavy_6_25MMDisp','ShearRect_6_25MMDisp'}
                    % If inputing your own Abaqus file
                    [x_temp,y_temp,z_temp,TRI_temp] = readInp([A{i} '_tet']);
                otherwise
                    [x_temp,y_temp,z_temp,TRI_temp] = genMesh(curdir,[A{i} '_tet'],mesh,params,'NoWrite',pres_disp{i},mesh_ref);
            end
        end
        x{i} = x_temp; y{i} = y_temp; z{i} = z_temp; TRI{i} = TRI_temp;
    end
else
    mesh = 'hex';
end

addpath(genpath('InputFiles'))

%% Main code
for ii = 1:length(A)
    % Change A{ii} to ['ShearWavy\' A{ii}] for folder organization
    switch mesh
        case 'hex'
            [X,Y,Z,ElNode,U] = runAbaqus(curdir,A{ii});
        case 'tet'
            [X,Y,Z,ElNode,U] = runAbaqus(curdir,[A{ii} '_tet'],x{ii},y{ii},z{ii},TRI{ii});
    end

    h = waitbar(0,'Progress: 0%');
    for idx = 1:length(ElNode)
        switch mesh
            case 'hex'
                for i = 1:8
                    n(i) = ElNode(idx,i);
                end
                du(1,:) = (U(n(5),:)-U(n(1),:)) + (U(n(8),:)-U(n(4),:)) + ...
                    (U(n(7),:)-U(n(3),:)) + (U(n(6),:)-U(n(2),:));
                du(2,:) = (U(n(1),:)-U(n(2),:)) + (U(n(4),:)-U(n(3),:)) + ...
                    (U(n(5),:)-U(n(6),:)) + (U(n(8),:)-U(n(7),:));
                du(3,:) = (U(n(1),:)-U(n(4),:)) + (U(n(2),:)-U(n(3),:)) + ...
                    (U(n(6),:)-U(n(7),:)) + (U(n(5),:)-U(n(8),:));
                du = du';
                dX(1) = X(n(5))-X(n(1)) + X(n(8))-X(n(4)) + X(n(7))-X(n(3)) + ...
                    X(n(6))-X(n(2));
                dX(2) = Y(n(1))-Y(n(2)) + Y(n(4))-Y(n(3)) + Y(n(5))-Y(n(6)) + ...
                    Y(n(8))-Y(n(7));
                dX(3) = Z(n(1))-Z(n(4)) + Z(n(2))-Z(n(3)) + Z(n(6))-Z(n(7)) + ...
                    Z(n(5))-Z(n(8));
                for i = 1:3
                    for j = 1:3
                        F{idx}(i,j) = du(i,j)/dX(j);
                    end
                end
                F{idx}=F{idx}+eye(3);
                F_pos{idx} = [mean(X(n(:))),...
                    mean(Y(n(:))),mean(Z(n(:)))];
                U_el{idx} = [mean(U(n(:),1)),...
                    mean(U(n(:),2)),mean(U(n(:),3))];
            case 'tet'
                for i = 1:size(ElNode,2)
                    n(i) = ElNode(idx,i);
                end
                dudN = [U(n(1),1)-U(n(4),1) U(n(2),1)-U(n(4),1) U(n(3),1)-U(n(4),1);
                    U(n(1),2)-U(n(4),2) U(n(2),2)-U(n(4),2) U(n(3),2)-U(n(4),2);
                    U(n(1),3)-U(n(4),3) U(n(2),3)-U(n(4),3) U(n(3),3)-U(n(4),3)];
                 dXdN = [X(n(1))-X(n(4)) X(n(2))-X(n(4)) X(n(3))-X(n(4));
                     Y(n(1))-Y(n(4)) Y(n(2))-Y(n(4)) Y(n(3))-Y(n(4));
                     Z(n(1))-Z(n(4)) Z(n(2))-Z(n(4)) Z(n(3))-Z(n(4))];
                 dudX = dudN*pinv(dXdN);
                 F{idx}=dudX+eye(3);
                 F_pos{idx} = [mean(X(n(:))),...
                    mean(Y(n(:))),mean(Z(n(:)))];
                 U_el{idx} = [mean(U(n(:),1)),...
                    mean(U(n(:),2)),mean(U(n(:),3))];
        end
        waitbar(idx/length(ElNode),h,...
            ['Progress for run #',num2str(ii),'/',num2str(length(A)),': ',num2str(floor(100*idx/length(ElNode))),'%'])
    end
    close(h)

    mod_met = min(cell2mat(F_pos'));
    for i = 1:3
        for idx = 1:length(ElNode)
            F_pos_test{i}(idx) = F_pos{idx}(i);
        end
        for j = 1:3
            F_t{ii}{i,j} = zeros(34,25,25);
            for idx = 1:length(ElNode)
                F_test{i,j}(idx) = F{idx}(i,j);
            end
        end
        for idx = 1:length(ElNode)
            U_test{i}(idx) = U_el{idx}(i);
        end
    end

    F_t{ii} = F_test;
    U_t{ii} = U_test;
    X_ref = F_pos_test;
end

if saveset
    switch mesh
        case 'hex'
            switch curdir
                case {'22-0201-5MM_Holes_SS','22-0301-NoHoles_SS','22-0325-Uniaxial',}
                    if ~exist(['Simulations\' curdir], 'dir')
                        mkdir(['Simulations\' curdir])
                    end
                    save(['Simulations\' curdir '\MRI-3Ddefs_SimpleShear_' curdir '.mat'], 'F_t', 'U_t');
                    save(['Simulations\' curdir '\refpositions.mat'],'X_ref');
                otherwise
                    if ~exist(['Simulations\' curdir], 'dir')
                        mkdir(['Simulations\' curdir])
                    end
                    save(['Simulations\' curdir '\MRI-3Ddefs_Iosipescu_' curdir '.mat'], 'F_t', 'U_t','lbl' );
                    save(['Simulations\' curdir '\refpositions.mat'],'X_ref');
            end
        case 'tet'
            switch curdir
                case {'22-0201-5MM_Holes_SS','22-0301-NoHoles_SS','22-0325-Uniaxial',...
                        '22-1212-Wavy_SS','22-1215-Wavy_Sweep','23-0119-Wavy_Prd_Sweep',...
                        '23-0202-Wavy-Thick_Amp_Sweep','23-0202-Wavy-Thick_Prd_Sweep',...
                        '23-0217-Wavy-Thick_Width_Sweep','23-0405-SingleCompression'}
                    if ~exist(['Simulations_tet\' curdir], 'dir')
                        mkdir(['Simulations_tet\' curdir])
                    end
                    save(['Simulations_tet\' curdir '\MRI-3Ddefs_SimpleShear_' curdir '.mat'], 'F_t', 'U_t');
                    save(['Simulations_tet\' curdir '\refpositions.mat'],'X_ref');
                otherwise
                    if ~exist(['Simulations_tet\' curdir], 'dir')
                        mkdir(['Simulations_tet\' curdir])
                    end
                    save(['Simulations_tet\' curdir '\MRI-3Ddefs_SimpleShear_' curdir '.mat'], 'F_t', 'U_t');
                    save(['Simulations_tet\' curdir '\refpositions.mat'],'X_ref');
                    save(['Simulations_tet\' curdir '\MRI-3Ddefs_Iosipescu_' curdir '.mat'], 'F_t', 'U_t','lbl' );
                    save(['Simulations_tet\' curdir '\refpositions.mat'],'X_ref');
            end
    end
end