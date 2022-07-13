function mshowspm(hdr, data)
%
% mshowspm(hdr,data)
% this function displays the SPM on top of whatever figure is
% being currently displayed. - Multiple slices version
%

global MAPMAX

	sz = size(data);
	matrix_size = round(sqrt(hdr.zdim));
	
	xdim = matrix_size * hdr.xdim;
	ydim = matrix_size * hdr.ydim;
	

%	for i=1:sz(1)
%		if data(i,3) == slice_number
%			sl = [sl; data(i,:)];
%		end
%	end	

	hold on ;

     spm_scale = MAPMAX/max(data(:,4))  ;
     
     
     
     for i=1:sz(1)
        
        row = ceil( data(i,3) / matrix_size);
        col = rem( data(i,3) , matrix_size) ;
        
        if col==0 
           col=matrix_size;
        end
        
        
        image(...
           data(i,1) + hdr.xdim*(row-1), ...
           data(i,2) + hdr.ydim*(col-1), ...
           data(i,4)*spm_scale + MAPMAX );
	end


	hold off

return













