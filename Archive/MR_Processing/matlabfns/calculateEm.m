function Em = calculateEm(Eij)
%Helper function for calculating von mises strain Em
%Input: gradient of displacement, E_ij{1:3,1:3}(Isize)
%Output: Von mises strain, Em(Isize)

Em = zeros(size(Eij{1,1}));
for i=3:-1:1
    for j=3:-1:1
        Em = Em + sqrt(2/3*(Eij{i,j}.*Eij{i,j}));
    end
end

end