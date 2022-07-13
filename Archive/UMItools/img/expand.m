function result=expand(support);
% function result=expand(support);
% expands the region by 1 pixel around
[M,N]=size(support);
result=1-shrink(1-support);


