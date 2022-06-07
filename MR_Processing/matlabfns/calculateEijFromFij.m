function Eij = calculateEijFromFij(Fij)
%Helper function for calculating Lagrangian strain tensor Eij
%Input: gradient of displacement, u_i,j{1:3,1:3}(Isize)
%Output: Lagrangian strain tensor, Eij{1:3,1:3}(Isize)

Iij = eye(3,3);

for i=3:-1:1
    for j=3:-1:1
        Eij{i,j} = 0.5*(Fij{j,i}.*Fij{i,j}-Iij(i,j));
    end
end

end