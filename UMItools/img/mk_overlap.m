function D = mk_overlap(D,inds,inds2)
msk1 = zeros(size(D));
msk1(inds)=1;
msk2 = zeros(size(D));
msk2(inds2)=1;
msk = msk1 .* msk2;
inds3 = find(msk);

if ~isempty(inds3)
    D(inds3) =  800;

end
return