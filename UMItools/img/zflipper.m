function zflipper(root, keepName)
% function zflipper(root, keepName)

if nargin<2
    keepName=0;
end

[d h] = read_img(root);
out = zeros(size(d));

for t=1:h.tdim
    
    d2 = d(t,:);
    d2 = reshape(d2,h.xdim, h.ydim, h.zdim);
    out2 = zeros(size(d2));
    
    for z=1:h.zdim
        out2(:,:,z) = d2(:,:,h.zdim-z+1);
    end
    
    out(t,:) = out2(:);
end

if keepName
    write_img([root '.img'], out, h);
else
    write_img(['r' root '.img'], out, h);
end

return
