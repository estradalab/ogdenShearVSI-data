function A = divv(mat,dim)
    if dim == 1
        A = mat(2:end,:,:)./mat(1:(end-1),:,:);
    elseif dim == 2
        A = mat(:,2:end,:)./mat(:,1:(end-1),:);
    elseif dim ==3
        A = mat(:,:,2:end)./mat(:,:,1:(end-1));
    end
end