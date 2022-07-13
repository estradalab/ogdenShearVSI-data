function [L Q V] = derivs2strain(M)

% FUNCTION derivs2strain converts a volume of displacement derivatives to a
%   volume of Lagrange strain tensors and a volume of equivalent strain
% 
%   Specifically written for use with MRI phase data (MR-u)
% 
%   INPUT: Matrix M, a volume of displacement derivatives
%           where M(i,j,k,:,:) = |du/dx du/dy du/dz|
%                                |dv/dx dv/dy dv/dz|
%                                |dw/dx dw/dy dw/dz|
%           and i,j,k are the indices corresponding to the x, y, z location
%               in the volume, respectively
%   OUTPUT: Matrix L, a volume of Lagrange strain tensors
%           where L(i,j,k,:,:) = |E11 E12 E13|
%                                |E21 E22 E23|
%                                |E31 E32 E33|
%       Matrix Q, a volume of equivalent strains where Q(i,j,k) is the
%           equivalent strain at position (i,j,k)
%       Matrix V, a volume of volumetric strains
% 
% Written by Callan Luetkemeyer
% University of Michigan
% 1/31/2019

[rows,cols,slices,tsiz1,tsiz2]=size(M);
L=zeros(rows,cols,slices,tsiz1,tsiz2);

for k=1:slices
    for j=1:cols
        for i=1:rows
            dUdX=squeeze(M(i,j,k,:,:));
            F=dUdX+eye(3);
            J=det(F);
            
            C=F'*F;
            
            Lprime=1/2*(C-eye(3));
            
            L(i,j,k,:,:)=Lprime;
            Edev=Lprime-1/3*trace(Lprime)*eye(3);
            Q(i,j,k)=sqrt(2/3*sum(dot(Edev,Edev)));
            V(i,j,k)=J-1;
        end
    end
end

end