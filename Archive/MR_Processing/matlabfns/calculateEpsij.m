function epsij = calculateEpsij(Gjui)
%Helper function for calculating Small strain tensor eps_ij
%Input: gradient of displacement, u_i,j{1:3,1:3}(Isize)
%Output: Small strain tensor, epsij{1:3,1:3}(Isize)

for i=3:-1:1
    for j=3:-1:1
        epsij{i,j} = 0.5*(Gjui{i,j}+Gjui{j,i});
    end
end

end