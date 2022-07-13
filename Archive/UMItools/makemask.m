function mask=makemask(in)
% this function generates a mask
% in is either a structure, or a 2D or 3D image
% if in is a structure, then a mask field will be added to it
% if in is a 2d or 3d matrix, then the output is a mask, on data that have
% been smoothed
voxsize=0.4; %blurring voxel is 0.4mm in size

if isstruct(in);
    if isfield(in,'image'); % this is a structure returned by varianms
        if ndims(in.image)==2;
            temp=blur(in.image);
            
        elseif ndims(in.image)>2;
            temp=blur3d(in,'vox',voxsize*[1 1 1]);
        end
        [signallevel, noiselevel ]=estimate_noiselevel(temp.image);
        mask=abs(temp.image)>5*noiselevel;
    end
elseif isnumeric(in);
    display('The input to makemask needs to be a structure made by varianms')
end


