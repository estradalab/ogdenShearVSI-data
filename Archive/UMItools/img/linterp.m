function f = interp(a,repfact)
% usage .. interp(a,f);
% increases matrix size of a by bilinear interpolation 
% by a factor of f (integer)

[szx, szy]  = size(a);
szxo = szx*repfact;
szyo = szy*repfact;
ahat = [a(2:szx,:)' zeros(szy,1)]';
% row space
aa = zeros(szxo,szy);
for j=1:repfact
  aa(j:repfact:szxo,:) = a.*((repfact-j+1)/repfact) + ahat.*((j-1)/repfact);
end
% column space
aahat = [aa(:,2:szy) zeros(szxo,1)];
f = zeros(szxo,szyo);
for j=1:repfact
  f(:,j:repfact:szyo) = aa.*((repfact-j+1)/repfact) + aahat.*((j-1)/repfact);
end
