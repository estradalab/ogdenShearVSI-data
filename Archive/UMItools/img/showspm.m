function showspm(data,slice_number)
%
% showspm(data, slice_number)
% this function displays the SPM on top of whatever figure is
% being currently displayed.
%
global MAPMAX 
global hdr

	sz = size(data);
   whos data
   
   MAPMAX;
	spm_scale = MAPMAX /max(data(:,4))  
	hold on ;
   
   	% Draw each one of the activated pixels at a time on top of the image...
      
  [i j] = ind2sub(sz , find(data(:,3)==slice_number));    
  sl = data(i,:);
  whos sl
    
  for i=1:size(sl,1)
         image(sl(i,1),sl(i,2), sl(i,4)*spm_scale + MAPMAX );
	end	

  
	hold off


return





function unused()



% this function would extract the SPM from the SPMt.mat file
% and convert them to pixels according to the header information.

global slices
global slice_number

	[file path] =uigetfile('*.mat', 'Select the SPM mat file');
	cd(path)
	load XYZ
	load (file)
	
	hdr = read_hdr('fmri001.hdr')

	
	PIX_mm = hdr.xsize;
	SLI_mm = hdr.zsize;
	X_org = hdr.xdim/2;
	Y_org = hdr.ydim/2;

	x_mm = XYZ(1,:)';
	y_mm = XYZ(2,:)';
	z_mm = XYZ(3,:)';

	x_pix = X_org - x_mm/PIX_mm;
	y_pix = Y_org - y_mm/PIX_mm;
	z_pix = z_mm/SLI_mm;
	
	data = [x_pix y_pix z_pix SPMt']

	oldmap = colormap;
	colormap = hot;

	sz = size(data)
	for i=1:sz(1)
		if (z_pix(i) == slice_number)
			imagesc(data(i,1), data(i,2), data(i,4))
		end
	end
	colormap = oldmap;

return












