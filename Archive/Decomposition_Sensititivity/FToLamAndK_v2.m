function [lambda, k] = FToLamAndK_v2(F)

C = F'*F;
[n, lam2] = eig(C);
Q = n;
%Note that n(:,1) is the first eigenvector, n(:,2) is the second, etc.

lam = sqrt(lam2);
lam = sort([lam(3,3) lam(2,2) lam(1,1)],'descend');

try

A = log(lam(1))/log(lam(2));
k_temp(1) = ((1.5+A) + sqrt((1.5+A)^2+2-2*A))/2;
k_temp(2) = ((1.5+A) - sqrt((1.5+A)^2+2-2*A))/2;
if sum(discretize(k_temp,[0,1]))==2
    k = [];
elseif discretize(k_temp(1),[0,1])==1
    k = k_temp(1);
elseif discretize(k_temp(2),[0,1])==1
    k = k_temp(2);
else
    k = [];
end
lambda = exp(log(lam(2))/(0.5-k+eps));

catch
    
    syms K L;

    [solK,solL] = vpasolve(L^(-K^2+3/2*K+1/2)==lam(1), L^(-K+1/2)==lam(2));%L^(-K+1/2)==lam(2)); (K^2-1/2*K-1)*log(L)==log(lam(3)) 

    k = eval(solK);
    lambda = eval(solL);
    
end

end

