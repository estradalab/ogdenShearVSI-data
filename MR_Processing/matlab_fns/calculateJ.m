function [J,Fbar] = calculateJ(Fij)

sz = size(Fij{1,1});
J = zeros(sz);
for a=1:sz(1)
    for b=1:sz(2)
        for c = 1:sz(3)
            F_ = [Fij{1,1}(a,b,c) Fij{1,2}(a,b,c) Fij{1,3}(a,b,c);...
                Fij{2,1}(a,b,c) Fij{2,2}(a,b,c) Fij{2,3}(a,b,c);...
                Fij{3,1}(a,b,c) Fij{3,2}(a,b,c) Fij{3,3}(a,b,c);];
            try
                J(a,b,c) = det(F_);
            catch
                J(a,b,c) = nan;
            end
        end
    end
end

for i=1:3
    for j=1:3
        Fbar{i,j} = Fij{i,j}.*J.^(-1/3);
    end
end

end