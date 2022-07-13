function dwcall(action)

% function dwcall(action)
%
% Luis Hernandez
% Last edit 1-9-97
%
% This function will hold all the different callback routines
% for the figure defined in dispwin.m
%
% The parameter 'action' will be the flag that indicates which 
% widget called the function and what it is that needs to be done.
%

global volume_data
global slice_number
global spm_overlay
global spm_data
global img_scale
global spm_scale
global imgname
global hdrname
global hdr
global multislice
global MAPMAX
global pix_sample
global AxisHandle 
global SPM_scale_factor


switch(action)
   
case 'load_*img'

   spm_overlay = 0;
   multislice = 1;
   pix_sample = 0;
   colormap(gray);

   s = size(colormap)
   MAPMAX =s(1);
   
   % Let user select filename ...
   [file path] = uigetfile('*.img','Select Analyze file');
   imgname = strcat(path,file);
   
   sz = size(imgname);
   hdrname = strcat(imgname(1:(sz(2) - 4)), '.hdr');
   hdr = read_hdr(hdrname);
   
   cd(path)
   
   slice_number = 1;
   volume_data= read_img2(hdr, imgname);
   
   disp ('read data successfully')
   
   mx = max(volume_data);
   img_scale = MAPMAX / max(max(mx));
   
   multisliceBox = findobj('Tag','msToggle');
   multislice = get(multisliceBox, 'Value')

   switch(multislice)
   case 0
      
      imdisp(volume_data(:,:,slice_number), img_scale);
      
   case 1
      mimdisp(hdr, img_scale, volume_data);
   end
      
   h= findobj('Tag','ImageFileName');
   set(h,'String',imgname);
   SliceMarker = findobj('Tag','EditSliceNum');
   set(SliceMarker, 'String', 1 );
   
   AxisHandle = gca
   
case 'multislice_toggle'
   %whos
   %[slice_number img_scale]

   
   multisliceBox = findobj('Tag','msToggle');
   multislice = get(multisliceBox, 'Value')
   
   if ~isempty(imgname)
      switch(multislice)
      case 0
         
         imdisp(volume_data(:,:,slice_number), img_scale);
         
         if spm_overlay == 1
            showspm(spm_data,slice_number);
         end
	
	   case 1
	 
         mimdisp(hdr, img_scale, volume_data);
         
         if spm_overlay == 1
            mshowspm(hdr, spm_data);
         end

	   end
	   
	   	 
      
     
   end
     
   
case 'load_spmt'
   %slices = spmimg('t');
   %slice_number = 1;
   %imdisp(slices(slice_number),img_scale);
   
case 'load_spmF'
   %slices = spmimg('F');
   %slice_number = 1;
  % imdisp(slices(slice_number), img_scale);
   
case 'save'
   %[file path] = uiputfile('*.img','Select Output File name');
   %imgname = strcat(path,file);
   %sz = size(imgname);
   %hdrname = strcat(imgname(1:(sz(2) - 4)), '.hdr');
   
   %hdr = make_hdr(slices);
   %write_hdr(hdrname, hdr);
   %write_img(slices, imgname);
    
case 'exit'
   
   clear slice;
   clear slice_number;
   main_window = findobj('Tag', 'MainDispWinTag')
   close (main_window);
   
case 'forward'
   if slice_number < hdr.zdim
      
      %cla
         slice_number = slice_number +1;
         
         imdisp(volume_data(:,:,slice_number), img_scale);
         
         SliceMarker = findobj('Tag','EditSliceNum');
         set(SliceMarker, 'String', num2str(slice_number) );
whos
         if spm_overlay == 1
            showspm(spm_data,slice_number);
         end

     end      
   %end
   
   
case 'changeslice'
   %cla
   SliceMarker = findobj('Tag','EditSliceNum');
   tmp = str2num(get(SliceMarker,'String'));
   if tmp <= hdr.zdim
      
      slice_number = tmp
      
      imdisp(volume_data(:,:,slice_number), img_scale);
      
      if spm_overlay == 1
         showspm(spm_data,slice_number);
      end
      
   end
   

case 'back'
   if slice_number > 1
      %cla
      slice_number = slice_number -1;
      
      imdisp(volume_data(:,:,slice_number), img_scale);
         
      SliceMarker = findobj('Tag','EditSliceNum');
      set(SliceMarker, 'String', num2str(slice_number) );
      
      if spm_overlay == 1
         showspm(spm_data,slice_number);
      end
      
   end
   
   
case 'timeseries'
   [x y] = ginput(2)
   xx  = [x(1) ; x(2) ; x(2) ; x(1); x(1) ];
   yy =  [y(1) ; y(1) ; y(2) ; y(2); y(1) ];
   line(xx,yy);

   z =slice_number;
   x = fix(x)
   y = fix(y)
   timeplot(x,y,z)


   
 case 'addspm'
    
   if spm_overlay == 0
      c = [gray;hot];  
      colormap(c);
      spm_overlay = 1;
      
   end
   
   %clear spm_data
   spm_data = extractspm2;
          
   switch(multislice)
   case 0
      showspm(spm_data, slice_number);

   case 1
      mshowspm(hdr, spm_data);
        
     end
 
    
   
case 'pixel_count'
      
   [x y] = ginput(2)
   x = fix(x)
   y = fix(y)
   
   line( x,y);
   
   midplane = [...
         x(1) y(1) 1;
         x(2) y(2) 1;
         x(2) y(2) 2]

      weighted_sums=zeros(2,5);
      
      if ~isempty(spm_data)
         sz = size(spm_data);
         bufsize = 100;
         
         i=1;
         while i< sz(1)
            if i+bufsize < sz(1)
               weighted_sums = weighted_sums + rightleft3( spm_data( i:i+bufsize-1, :),midplane);
               i=i+bufsize;
            else
               weighted_sums = weighted_sums + rightleft3( spm_data( i:sz(1), :), midplane);
               i=sz(1);
            end
            
            % recall that rightleft3 returns 
            % [x_right y_right z_right total_right_weight rightcount;
            %  x_left  y_left  z_left total_left_weight leftcount]
            
         end
         
      else
         errormesg('Need to add SPM first');
      end
      
   
   
      weighted_aves = zeros(2,3)
      for i=1:2
         weighted_aves(i,:) = weighted_sums(i,1:3) / weighted_sums(i,4);
      end
      
      disp('right and left centroids')
      weighted_aves
      
      disp('right and left counts')
      weighted_sums(:,5)
      plot3(weighted_aves(:,1),weighted_aves(:,2),weighted_aves(:,3),'g*')
      
     
      [fn,pn] = uiputfile('*.txt','Save Stats As ...');
      fn = strcat(pn,fn);
      
      fp = fopen(fn,'w');
      
      fprintf(fp,'voxel size:	\t%f\t%f\t%f\n\n',hdr.xsize, hdr.ysize, hdr.zsize);

      fprintf(fp,'Right voxels:	\t%d\n',weighted_sums(1,5));
      fprintf(fp,'Left  voxels:	\t%d\n',weighted_sums(2,5));
      fprintf(fp,'Right Centroid: \t%f\t%f\t%f\n',weighted_aves(1,:) );
      fprintf(fp,'Left  Centroid: \t%f\t%f\t%f\n',weighted_aves(2,:) );
   
      fclose (fp);
   
  
   case 'get_coords'
      
      handle = findobj('Tag','CoordsButton');
      pix_sample = get(handle,'Value')
      
      matrix = ceil(sqrt(hdr.zdim));
   
      while pix_sample==1
         
         disp('Looping...')
         xy = ginput(1)
         CoordinatesDisp = findobj('Tag','pixel_coordinates');

         if (multislice)
            xyz(1) = rem(xy(1), hdr.xdim);
            xyz(2) = rem(xy(2), hdr.ydim);
            xyz(3) = floor(xy(1)/hdr.xdim) *matrix + floor(xy(2)/hdr.ydim) +1;
            val = volume_data(xyz(1),xyz(2),xyz(3));
            set(CoordinatesDisp, 'String', num2str([floor(xyz) val]) );
          
         else
            val = volume_data(xy(1),xy(2),slice_number) 
            set(CoordinatesDisp, 'String', num2str([floor(xy) slice_number val]) );
            xyz = [xy slice_number];
         end
         
         if spm_overlay
            xyz = floor(xyz);
            a = spm_data(find(spm_data(:,1)== xyz(1)),:);
            b = a(find(a(:,2)== xyz(2)),:);
            c = b(find(b(:,3)== xyz(3)),:)
            c(:,4) = c(:,4)* SPM_scale_factor
            
            if ~isempty(c)
               set(CoordinatesDisp, 'String', strcat(num2str([c]),'(SPM)') );
            end
               
         end
         
         handle = findobj('Tag','CoordsButton');
         pix_sample = get(handle,'Value')

         
      end
      

         
      
  case 'reset'
     cla
     clear img_scale slice slice_number spm_overlay spm_data spm_scale
     colormap(gray);
     
     
  end

return















