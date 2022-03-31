function defgrad = convDennisStoJonS_lin(S_t)

volsz = size(S_t); %should be t, i, j, k
%     F11 = ones(szDG(2),szDG(3),szDG(4));
%     F22 = F11; F33 = F11;
%     F12 = zeros(szDG(2),szDG(3),szDG(4));
%     F13 = F12; F21 = F12; F23 = F12; F31 = F12; F32 = F12;
%
%% Match up the indices by looping
try
for t = volsz(1):-1:1
    for i = volsz(2):-1:1
        for j = volsz(3):-1:1
                defgrad{t}(i,j) = S_t(t,i,j);
        end
    end

end

catch
    for i = volsz(1):-1:1
        for j = volsz(2):-1:1
                defgrad(i,j) = S_t(i,j);
        end
    end

end

end