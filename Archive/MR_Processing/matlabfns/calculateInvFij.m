function Finv_ij = calculateInvFij(Fij)
%Helper function for calculating Lagrangian strain tensor Eij
%Input: gradient of displacement, u_i,j{1:3,1:3}(Isize)
%Output: Lagrangian strain tensor, Eij{1:3,1:3}(Isize)

defgrad = convJonFtoCalF(Fij);
volsz = size(defgrad);

for i=volsz(1):-1:1
    for j=volsz(2):-1:1
        for k=volsz(3):-1:1
            %Eij{i,j} = 0.5*(Fij{j,i}.*Fij{i,j}-Iij(i,j));
            if ~isnan(defgrad{i,j,k}(1,1))
                Finv_{i,j,k} = inv(defgrad{i,j,k});
            else
                Finv_{i,j,k} = nan(3,3);
            end
        end
    end
end

Finv_ij = convCalFtoJonF(Finv_);

end