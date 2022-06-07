function outaxis=newaxis(oldaxis,numberofentries)

spacevec=diff(oldaxis);
if any(abs(spacevec/spacevec(1)-1)>1e-10);
	display('newaxis: wonky input axis has uneven spacing');
end
spacing=mean(spacevec);
FOV=spacevec(1)*numel(oldaxis);
Limit1=oldaxis(1)-spacing/2;      %this accounts for the voxel size
Limit2=oldaxis(end)+spacing/2;

if abs(abs(FOV)-abs(Limit1-Limit2))/abs(FOV)>1e-10;
	display('Error in axis calculation');
	return
end
newspacing=(Limit2-Limit1)/numberofentries;
outaxis=Limit1+[0.5:1:numberofentries-0.5]*newspacing;

end


