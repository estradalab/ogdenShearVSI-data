function f=pixtoim(pix_list,im_size,x_list,filter)
% Usage ... f=pixltoim(pix_list,im_size,x_list,filter)
%
% x_list=data (vector), pix_list=corresponding pixels(2 col)
% x_list and pix_list must be the same length
% x_list must be a vector and pix_list is a 2col vector
% filter is a binary vector which indicates inclusion/exclusion

if nargin<4, filter=ones([length(pix_list) 1]); end;

if nargin<3,
  x_list=ones(size(pix_list));
end;

if nargin<2,
  tmp=max(log2(pix_list));
  imsize=[2^(floor(tmp(1))) 2^(floor(tmp(2)))];
end;

f=zeros(im_size);

for m=1:size(x_list,1),
  f(pix_list(m,1),pix_list(m,2))=f(pix_list(m,1),pix_list(m,2))+ ...
                                   filter(m)*x_list(m);
end;

if nargout==0,
  show(f');
end;

