function RT = readair(AP,V,u)
% function RT = readair(AP,V,u)
%
% Read Rotations & Translations from AIR file
% FORMAT RT = readair(AP)
% AP  -  AIR parameter file (".air file")
% V   -  Verbosity
% u   -  Units ("mm", "pixels", or "ASmm")
% RT  -  [pitch roll yaw xshift yshift zshift]
%___________________________________________________________________________
%
% Rotations are in units of degrees.
% Translations are in real units of the images (usually millimeters).
%
% If V is 1, rotations and translations are printed nicely.
% If V is 2, rotations and translations are printed in a form that
% matches Anne Smith's registration programs.
%
% If a filename of '0' is specified, RT = [0 0 0 0 0 0]
%
% This uses the formulation of the AIR file found on the 
% Automated Image Registration Web page.
%
% If there are spaces in AP, multiple files are assumed.
%
% "mm" are the default units.
%
% @(#)readair.m	2.4 Tom Nichols, UPMC PET Facility 96/07/24


if (nargin < 2); V = 0; end
ConvUnits = 1; Units = 'mm';
if (nargin >= 3)
    if strcmp(lower(u),'pixels')
	ConvUnits = 0;
	Units = 'pixels';
    end
end

RT=[];
[Nm AP] = strtok(AP,setstr([9,10,13,32]));
Nm = Nm(Nm ~= 10);  % strtok doesn't recognize NL

while (~isempty(Nm))

    if (strcmp(Nm,'0'))
	RT = [RT;0 0 0 0 0 0];
    else
    
        %
        % Open file
        %
	fid     = fopen(Nm);
	if (fid <= 0);
	    RT = [];
	    return
	end

	%
	% Check that it's an AIR file	
	%	
	fseek(fid,0,'eof');
	len=ftell(fid);
	fseek(fid,0,'bof');	
	if len~=688 & len~=808 & len~=720
	    fprintf(2,'''%s'' is not an AIR file',Nm);
	    RT = [];
	    return
	end	    
	if len==720
	    AIRver = 2;
	else
	    AIRver = 1;
	end
	    
        %    
        % Read in all of the AIR file, though we don't use all of it &
        % build E, the Homogenious transformation matrix.
	%	
	if AIRver==2	
	    Eair    = fread(fid,[4 4],'double');
	    E       = [Eair];
	else	
	    Eair    = fread(fid,[3 4],'double');
	    % Old AIR used funny order and had no perspective 
	    E       = [Eair(:,2:4) Eair(:,1);0 0 0 1];
	end	    
	s_file  = setstr(fread(fid,128,'uchar')');
	tmp     = fread(fid,4,'int');
	s_bits  = tmp(1); sx_dim = tmp(2); sy_dim = tmp(3); sz_dim = tmp(4);
	tmp     = fread(fid,3,'double');
	sxSz = tmp(1); sySz = tmp(2); szSz = tmp(3);
	r_file  = setstr(fread(fid,128,'uchar')');
	tmp     = fread(fid,4,'int');
	r_bits = tmp(1); rx_dim = tmp(2); ry_dim = tmp(3); rz_dim = tmp(4);
	tmp     = fread(fid,3,'double');
	rxSz = tmp(1); rySz = tmp(2); rzSz = tmp(3);
	comment = setstr(fread(fid,128,'uchar')');
	tmp     = fread(fid,2,'ulong');
	s_hash = tmp(1); r_hash = tmp(2);
	tmp     = fread(fid,2,'ushort');
	s_volume = tmp(1); r_volume = tmp(2);
	reserved = fread(fid,116,'uchar');
	
	fclose(fid);

	%
        % E matrix components...
        %
        % E = Zr * Cr * T * R * P * Cs
	%
	% [x_r;y_r;z_r;1] = E * Zs * [x_s;y_s;z_s;1]
	%
	%
        %
	if AIRver == 2	
	    s_VxSz = min([sxSz,sySz,szSz]);
	    sx_zoom = sxSz/s_VxSz;
	    sy_zoom = sySz/s_VxSz;
	    sz_zoom = szSz/s_VxSz;
	    r_VxSz = min([rxSz,rySz,rzSz]);
	    rx_zoom = rxSz/r_VxSz;
	    ry_zoom = rySz/r_VxSz;
	    rz_zoom = rzSz/r_VxSz;
	    
	    Zs = diag([sx_zoom;sy_zoom;sz_zoom;1]);
	    Cs = [1 0 0 -(sx_dim-1)*sx_zoom/2; 0 1 0 -(sy_dim-1)*sy_zoom/2; 0 0 1 -(sz_dim-1)*sz_zoom/2;0 0 0 1];
	    P = [sxSz/rxSz 0 0 0; 0 sxSz/rxSz 0 0; 0 0 sxSz/rxSz 0; 0 0 0 1];
	    % We don't know R & T yet
	    R = zeros(4,4);
	    T = zeros(4,4);
	    Cr = [1 0 0 (rx_dim-1)*rx_zoom/2; 0 1 0 (ry_dim-1)*ry_zoom/2; 0 0 1 (rz_dim-1)*rz_zoom/2;0 0 0 1];
	    Zr = diag([1/rx_zoom;1/ry_zoom;1/rz_zoom;1]);
	    
	else
	    szoom=szSz/sxSz;
	    rzoom=rzSz/rxSz;
	    r_VxSz = rxSz;	    
	    
	    Zs = diag([1;1;szoom;1]);
	    Cs = [1 0 0 -(sx_dim-1)/2; 0 1 0 -(sy_dim-1)/2; 0 0 1 -(sz_dim-1)*szoom/2;0 0 0 1];
	    P=[sxSz/rxSz 0 0 0; 0 sxSz/rxSz 0 0; 0 0 sxSz/rxSz 0; 0 0 0 1];
	    % We don't know R & T yet
	    R = zeros(4,4);
	    T = zeros(4,4);
	    Cr = [1 0 0 (rx_dim-1)/2; 0 1 0 (ry_dim-1)/2; 0 0 1 (rz_dim-1)*rzoom/2;0 0 0 1];
	    Zr = diag([1;1;1/rzoom;1]);
	    
	end	    

	TandR = Cr\(Zr\(E/Cs));
    
        % E goes from Std to Resl, but you specify movements from Resl to Std
        % (e.g. w/ manualreslice), so we need minus sign here... (rotations 
	% don't have unambiguous direction, so there's no correction)
	T = -TandR(1:3,4);
	R = TandR(1:3,1:3);

        % Figure out rotations
	cosPitch = sqrt(R(2,1)^2 + R(2,2)^2);
	Pitch	 = atan2(R(2,3),cosPitch);
	if cosPitch == 0 
	    YawLesRoll = atan2(R(3,1),R(1,1));
	    Yaw = YawLesRoll/2*sign(YawLesRoll);
	    Roll = Yaw - YawLesRoll;
	else	
	    Roll = atan2(R(1,3)/cosPitch,R(3,3)/cosPitch);
	    Yaw	 = atan2(-R(2,1)/cosPitch,R(2,2)/cosPitch);
	end	    

	xshift	= T(1);
	yshift	= T(2);
	zshift	= T(3);
    
	% Add to RT, converting rotations to degrees, translations to real
	% units
	if (ConvUnits)
	    RT = [RT; [Pitch Roll Yaw]/(2*pi)*360 [xshift yshift zshift]*r_VxSz];
	else
	    RT = [RT; [Pitch Roll Yaw]/(2*pi)*360 [xshift yshift zshift]];
	end
	
	if (V==1)
	    % alignpettopet-style verbose output
	    rt = RT(size(RT,1),:);
	    n = sprintf('%s:\n',Nm);
	    a = sprintf('pitch=%f degrees roll=%f degrees yaw=%f degrees\n',rt(1:3));
	    b = sprintf('x_shift=%f %s y_shift=%f %s z_shift=%f %s\n', ...
		rt(4),Units,rt(5),Units,rt(6),Units);
	    disp([n a b])
	elseif (V==2)	    
	    % alignlinear-AnneSmith-mod-style verbose output
	    rt = RT(size(RT,1),:);
	    n = sprintf('%s:',Nm);
	    a = sprintf('%f,%f,%f,%f,%f,%f\n', ...
		    rt([3,1,2,4,5,6]));
	    disp([n a])
	end
    end
    
    [Nm AP] = strtok(AP,' 	');
    if ~isempty(Nm); Nm = Nm(Nm ~= 10); end
end


return



