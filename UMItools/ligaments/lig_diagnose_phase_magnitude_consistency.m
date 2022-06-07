function out=lig_diagnose_phase_magnitude_consistency(in)
% takes the three phase encoding data sets and goes through them to check
% in the region where there is signal whether or not the phase behaviour is
% credible

out=in;

ND=numel(in);  %number of data sets

for jj=1:ND
	
	M=in(jj).image; % (rd,ph,sl,encodepolarity)
	ph=M./abs(M);		%phase factor
	mm=mean(abs(M),4);	%mean amplitude
	pph=prod(ph,4);	%product of phase factors
	
	
	
	
	
	
	out.prodenc=mm.*pph;
	
	[sig, noise]=estimate_noiselevel(mm); %#ok<ASGLU> %get signal and noise level
	mask=abs(mm)>2*noise; %make a mask
	
	%find how many valid pixels there are in each slice
	si=size(mask);
	for sl=1:si(3)
		dummy=mask(:,:,sl);
		wt(sl)=mean(dummy(:));
	end
	
	%pick the slices with lots of signal
	wt=wt/max(wt); slind=find(wt>0.75);
	% shim those slices
	for sl=1:numel(slind)
		[outmatrix,rotmatrix{sl},aux(sl)]=gradshim(pph(:,:,slind(sl)));
		%figure; for jj=1:16; mesh(aux(jj).bowl+jj*5); hold on; end;  %test
	end
	test=cat(3,aux(:).bowl)/2;
	bowl(:,:,jj)=mean(test,3);
end

superbowl=mean(bowl,3);  
superphase=exp(-1i*superbowl);  %common phase factor for each pixel in a plane

%rotate the data for visualization
siM=size(M);
for jj=1:ND
	for sli=1:siM(3);
		out(jj).image(:,:,sli,1)=in(jj).image(:,:,sli,1).*superphase;
		out(jj).image(:,:,sli,2)=in(jj).image(:,:,sli,2).*superphase;
	end
end
	

end

