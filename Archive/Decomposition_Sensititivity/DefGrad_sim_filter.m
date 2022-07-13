function [F_t, U_t] = DefGrad_sim_filter(F_t,U_t)
for t = 1:2
    for i = 1:3
        for j = 1:3
            mask = F_t{t}{i,j};
            mask(~isnan(mask)==1) = 1;
            F_t{t}{i,j} = ndnanfilter(F_t{t}{i,j},'hamming',[2 2 2]).*mask;
        end
        U_t{t}{i} = ndnanfilter(U_t{t}{i},'hamming',[2 2 2]).*mask;
    end
end
end