function [E_t] = F_t_to_E_t(F_t)
F_temp = convJonFtoDennisF(F_t);
tlen = size(F_t,2);
volsz = size(F_t{1,tlen}{1,1});
for t = tlen:-1:1
    for i = volsz(1):-1:1
        for j = volsz(2):-1:1
            for k = volsz(3):-1:1
                F_calc = F_temp{t,i,j,k};
                E_calc{t,i,j,k} = 0.5*(transpose(F_calc)*F_calc-eye(3));
            end
        end
    end
end
E_t = convDennisFtoJonF(E_calc);
end