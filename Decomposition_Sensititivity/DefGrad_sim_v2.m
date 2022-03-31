clc; clear;

%% Sample with No Holes
% curdir = '22-0131-OgdenSS_Simulation';
% addpath(genpath(curdir))
% cd(curdir)
% xlRange = 'B2:G23661';
% xlRange2 = 'B2:I21251';
% A{1} = xlsread('OgdenNeoHookeRubber.xlsx','Ogden',xlRange);
% ElNode = xlsread('OgdenNeoHookeRubber.xlsx','ElNodes',xlRange2);
% holes = false;

%% Sample with Holes
% curdir = '22-0201-OgdenSS_Simulation_Holes';
% addpath(genpath(curdir))
% cd(curdir)
% xlRange = 'B2:G32621';
% xlRange2 = 'B2:I29275';
% A{1} = xlsread('OgdenNeoHookeRubber_2_5MM.xlsx','Ogden',xlRange);
% A{2} = xlsread('OgdenNeoHookeRubber_5MM.xlsx','Ogden',xlRange);
% ElNode = xlsread('OgdenNeoHookeRubber_2_5MM.xlsx','ElNodes',xlRange2);
% holes = true;

%% Sample with Holes (HIRES)
% curdir = '22-0209-OgdenSS_Simulation_Holes_HIRES';
% addpath(genpath(curdir))
% cd(curdir)
% xlRange = 'B2:G49316';
% xlRange2 = 'B2:I44847';
% A{1} = xlsread('OgdenNeoHookeRubber_2_5MM.xlsx','Ogden',xlRange);
% A{2} = xlsread('OgdenNeoHookeRubber_5MM.xlsx','Ogden',xlRange);
% ElNode = xlsread('OgdenNeoHookeRubber_2_5MM.xlsx','ElNodes',xlRange2);
% holes = true;

%% Sample with No Holes (3 Load Cells)
curdir = '22-0301-OgdenSS_Simulation_NoHoles_Corr';
addpath(genpath(curdir))
cd(curdir)
xlRange = 'B2:G23661';
xlRange2 = 'B2:I21251';
A{1} = xlsread('OgdenNeoHookeRubber_2_5MM.xlsx','Ogden',xlRange);
A{2} = xlsread('OgdenNeoHookeRubber_5MM.xlsx','Ogden',xlRange);
A{3} = xlsread('OgdenNeoHookeRubber_7MM.xlsx','Ogden',xlRange);
ElNode = xlsread('OgdenNeoHookeRubber_2_5MM.xlsx','ElNodes',xlRange2);
holes = false;

for ii = 1:length(A)
    X = A{ii}(:,1); Y = A{ii}(:,2); Z = A{ii}(:,3); % Coordinates
    U = A{ii}(:,4:6); % Displacements

    h = waitbar(0,'Progress: 0%');
    for idx = 1:length(ElNode)
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
                F{idx}(i,j)= du(i,j)/dX(j);
            end
        end
        F{idx}=F{idx}+eye(3);
        F_pos{idx} = [mean(X(n(:))),...
            mean(Y(n(:))),mean(Z(n(:)))];
        U_el{idx} = [mean(U(n(:),1)),...
            mean(U(n(:),2)),mean(U(n(:),3))];
        waitbar(idx/length(ElNode),h,...
            ['Progress: ',num2str(floor(100*idx/length(ElNode))),'%'])
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

    F11_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',F_test{1,1}');
    F12_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',F_test{1,2}');
    F13_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',F_test{1,3}');
    F21_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',F_test{2,1}');
    F22_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',F_test{2,2}');
    F23_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',F_test{2,3}');
    F31_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',F_test{3,1}');
    F32_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',F_test{3,2}');
    F33_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',F_test{3,3}');

    U1_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',U_test{1}');
    U2_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',U_test{2}');
    U3_int = scatteredInterpolant(F_pos_test{1}',F_pos_test{2}',F_pos_test{3}',U_test{3}');

    F_int = {F11_int F12_int F13_int; F21_int F22_int F23_int; F31_int F32_int F33_int};
    U_int = {U1_int,U2_int,U3_int};

    [s_x,s_y,s_z] = ndgrid(min(F_pos_test{1}):(max(F_pos_test{1})-min(F_pos_test{1}))/33:max(F_pos_test{1}),...
        min(F_pos_test{2}):(max(F_pos_test{2})-min(F_pos_test{2}))/24:max(F_pos_test{2}),...
        min(F_pos_test{3}):(max(F_pos_test{3})-min(F_pos_test{3}))/24:max(F_pos_test{3}));

    sz = size(F_t{1}{3,1});
    for i = 1:3
        for j = 1:3
            F_int_temp = F_int{i,j};
            F_t_temp = F_int_temp(s_x,s_y,s_z);
            if holes
                F_t_temp(:,10:15,[7:9 16:18]) = NaN;
                F_t_temp(:,12:13,[6 10 15 19]) = NaN;
                F_t_temp(:,[11 14],[10 15]) = NaN;
            end
            F_t{ii}{i,j} = F_t_temp;
        end
        U_int_temp = U_int{i};
        U_t_temp = U_int_temp(s_x,s_y,s_z);
        if holes
            U_t_temp(:,10:15,[7:9 16:18]) = NaN;
            U_t_temp(:,12:13,[6 10 15 19]) = NaN;
            U_t_temp(:,[11 14],[10 15]) = NaN;
        end
        U_t{ii}{i} = U_t_temp;
    end
end


cd ..