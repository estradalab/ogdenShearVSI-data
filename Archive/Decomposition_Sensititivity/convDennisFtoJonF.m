function defgrad = convDennisFtoJonF(F_t)

volsz = size(F_t); %should be t, i, j, k
try
%     F11 = ones(szDG(2),szDG(3),szDG(4));
%     F22 = F11; F33 = F11;
%     F12 = zeros(szDG(2),szDG(3),szDG(4));
%     F13 = F12; F21 = F12; F23 = F12; F31 = F12; F32 = F12;
%
%% Match up the indices by looping
for t = volsz(1):-1:1
    for i = volsz(2):-1:1
        for j = volsz(3):-1:1
            for k = volsz(4):-1:1
                defgrad{t}{1,1}(i,j,k) = F_t{t,i,j,k}(1,1);
                defgrad{t}{1,2}(i,j,k) = F_t{t,i,j,k}(1,2);
                defgrad{t}{1,3}(i,j,k) = F_t{t,i,j,k}(1,3);
                defgrad{t}{2,1}(i,j,k) = F_t{t,i,j,k}(2,1);
                defgrad{t}{2,2}(i,j,k) = F_t{t,i,j,k}(2,2);
                defgrad{t}{2,3}(i,j,k) = F_t{t,i,j,k}(2,3);
                defgrad{t}{3,1}(i,j,k) = F_t{t,i,j,k}(3,1);
                defgrad{t}{3,2}(i,j,k) = F_t{t,i,j,k}(3,2);
                defgrad{t}{3,3}(i,j,k) = F_t{t,i,j,k}(3,3);
            end
        end
    end

end

catch
    for i = volsz(2):-1:1
        for j = volsz(3):-1:1
            for k = volsz(4):-1:1
                defgrad{1,1}(i,j,k) = F_t{i,j,k}(1,1);
                defgrad{1,2}(i,j,k) = F_t{i,j,k}(1,2);
                defgrad{1,3}(i,j,k) = F_t{i,j,k}(1,3);
                defgrad{2,1}(i,j,k) = F_t{i,j,k}(2,1);
                defgrad{2,2}(i,j,k) = F_t{i,j,k}(2,2);
                defgrad{2,3}(i,j,k) = F_t{i,j,k}(2,3);
                defgrad{3,1}(i,j,k) = F_t{i,j,k}(3,1);
                defgrad{3,2}(i,j,k) = F_t{i,j,k}(3,2);
                defgrad{3,3}(i,j,k) = F_t{i,j,k}(3,3);
            end
        end
    end

end

end