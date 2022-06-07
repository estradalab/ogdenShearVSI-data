function [data h ] = fdf2nii(pathname, doJPG)

% function [data h ] = fdf2nii(pathname, doJPG)

% This script reads in Varian FDF files (2D and 3D)
% Use: type ReadVarFDF at the matlab command prompt
% Currently 3D data is displayed through the middle of the 3D dimension
% rank1 matrix1 storage1 bits1 type1
% rank is the dimension of the file
% matrix is the number of data points in each directions
%
% storage (float, double, int), etc (string)-so far all are float32 type
% bits is an integer
%
% written by Lana Kaiser, Varian, Dec. 2009
%
% Luis's addition:  read the whole volume and create a nifti header and
% data block
%


fnames = dir([pathname  'slice*']);
h = define_nii_hdr;

filename = fnames(1).name;
if filename ~=0
	%Open the input file
	filename=[pathname filename];
	[fid,msg]=fopen(filename,'r');
	if fid<0
		% There is an error
		disp(['File Open Failed: ' filename])
		return
	else
		disp(['reading ....' filename])
	end
end

%Start reading the header of the file
counter_line = 0;
line = '  ';

while  (counter_line < 100)


	line = fgetl(fid);
	%Get the important parameters: rank, matrix, format:

	if strmatch('float  rank = ', line) %determine 2D or 3D
		rank1_index=isstrprop(line, 'digit');
		rank1 = str2num(line(find(rank1_index>0)));

	end


	if strmatch('float  bits = ', line)
		eqsign=findstr(line,'=');
		scolon=findstr(line,';');

		bits1 = str2num(strtrim(line(eqsign+1:scolon-1)));

	end

	if strmatch('float  roi[] = ', line)
		comma=findstr(line,',');
		bracket1=findstr(line, '{');
		bracket2=findstr(line, '}');
		xfov = str2num(line(bracket1+1:comma(1)-1));
		yfov = str2num(line(comma(1)+1:comma(2)-1));
		zfov = str2num(line(comma(2)+1:bracket2-1));
	end

	if strmatch('float  matrix[] = ', line) %size of the total matrix
		if rank1==2; %this is 2D case

			matrix1_index=isstrprop(line, 'digit');
			comma=findstr(line,',');
			totaldigits = (find(matrix1_index>0));
			dim(1) = str2num(line(totaldigits(1):comma-1));
			dim(2) = str2num(line(comma+2:totaldigits(end)));

		elseif rank1==3;

			matrix1_index=isstrprop(line, 'digit');
			comma=findstr(line,',');
			totaldigits = (find(matrix1_index>0));

			dim(1) = str2num(line(totaldigits(1):comma(1)-1));
			dim(2) = str2num(line(comma(1)+2:comma(2)-1));
			dim(3) = str2num(line(comma(2)+2:totaldigits(end)))

		end
	end

	if strmatch('int    slices = ', line) %size of the total matrix
		pos = isstrprop(line, 'digit');
		dim(3) = str2num(line(pos));

	end

	if strmatch('int    bigendian', line)
		[token, rem] = strtok(line,'int   bigendian = ');
		endian = 'l'; % PC format (linux)
		if str2num(token)==1;
			endian = 'b';
		end

	end

	if strmatch('float  array_dim ', line)
		pos = find(line == '=');
		dim(4) = str2num(line(pos+1:end));
	end

	if strmatch('float  gap ', line)
		pos = find(line == '=');
		gap = str2num(line(pos+1:end));
	end


	counter_line=counter_line+1;

end

fclose(fid)

% if this is the first slice we're reading, create a header and
% allocate space now.

h.dim(1) = 4;
h.dim(2) = dim(1);
h.dim(3) = dim(2);
h.dim(4) = dim(3);
h.dim(5) = dim(4);

h.pixdim(2) = xfov / dim(1);
h.pixdim(3) = yfov / dim(2);
h.pixdim(4) = zfov + gap ;

h.bitspix = bits1;
h.datatype = 16;

data = zeros(h.dim(2), h.dim(3), h.dim(4), h.dim(5) );

%We are done with figuring out the file size info, now move on to read the
%binary portion
%%

for f=1:dim(4)  % images (frames in the movie)

	for s=1:dim(3)  % slices in each frame

		filename = sprintf('%sslice%03dimage%03decho001.fdf',pathname, s, f);

		[fid,msg] = fopen(filename,'r');

		if fid<0
			% There is an error
			disp(['File Open Failed: ' filename])
			return
		else
			%disp(['reading ....' filename])
		end

		if rank1==2; %this is 2D case
			%disp('reading 2D...')
			status = fseek(fid, -dim(1)*dim(2)*bits1/8, 'eof');

			image_space=fread(fid,[dim(1), dim(2)],'float32',endian);

			data(:, :, s, f) = image_space;

			imagesc(image_space');
			axis image
			set(gca,'DataAspectRatio', [ xfov yfov 1]);
			colormap(gray)
			if doJPG
				print(gcf,'-djpeg',['image_' num2str(f,'%05d')]);
			end

		else
			%if rank1==3;
			%this is 3D case
			disp('reading 3D...')
			status = fseek(fid, -dim(1)*dim(2)*dim(3)*bits1/8, 'eof');


			image_space=fread(fid,[dim(1)*dim(2), dim(3)],'float32',endian);

			data (:,:,:,f) = reshape(image_space,[dim(1) dim(2) dim(3)]);

		end

		if status<0
			disp('Error reading file')
		end

		fclose(fid);


	end
end

data2 = (reshape(data,  h.dim(2)*h.dim(3)*h.dim(4), h.dim(5)))';

write_nii('tmp.nii', data2, h, 0);


return
