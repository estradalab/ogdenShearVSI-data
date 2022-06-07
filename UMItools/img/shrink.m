function result=shrink(support);
% function result=shrink(support);
% shrinks the region by 1 pixel around
[M,N]=size(support);
result=support;
for ii=1:M
for jj=1:N
if support(ii,jj)==0
result(max(ii-1,1):min(ii+1,M),max(jj-1,1):min(jj+1,N))=0;
end;
end;
end;

