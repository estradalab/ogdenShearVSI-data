function [k1,k2,k3] = FToK2AndK3(F)

C = F'*F;
[n, lam2] = eig(C);
Q = n;
%Note that n(:,1) is the first eigenvector, n(:,2) is the second, etc.

lam = sqrt(lam2);
lam = sort([lam(3,3) lam(2,2) lam(1,1)],'descend');

k1 = log(lam(1)*lam(2)*lam(3));
for i = 1:3
    g(i) = log(lam(i))-k1/3;
end
k2 = sqrt(g(1)^2+g(2)^2+g(3)^2);
k3 = 3*sqrt(6)*g(1)*g(2)*g(3)/k2^3;

end

