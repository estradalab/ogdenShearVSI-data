function Eij = calculateEij(Gjui)
%Helper function for calculating Lagrangian strain tensor Eij
%Input: gradient of displacement, u_i,j{1:3,1:3}(Isize)
%Output: Lagrangian strain tensor, Eij{1:3,1:3}(Isize)

for i=3:-1:1
    for j=3:-1:1
        Eij{i,j} = 0.5*(Gjui{i,j}+Gjui{j,i}...
            +Gjui{1,i}.*Gjui{1,j}+Gjui{2,i}.*Gjui{2,j}+Gjui{3,i}.*Gjui{3,j});
    end
end

end