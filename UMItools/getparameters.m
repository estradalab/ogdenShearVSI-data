function out=getparameters(procparfile)
%function out=getparameters(procparfile)
    parlist={...
        'seqcon';               %compressed/standard loop control
        'ident';                %animal identity
        'studyid_';              %study id
        'te';                   %echo time
        'te2';                  %echo spacing with multi echo
        'te3';                  %echo spacing with multi echo?
		'te0';                  %pilot image after the 90 in UMS se sequences
		'te00';                 %separation between s-echo and subsequent multiecho
        'esp';                  %echospacing in FSE and the like
        'tr';                   %repeat time
        'ti';                   %inversion time
        'at';                   %acquisition time
        'gain';                 %receiver gain
        'flip1';                %excitation pulse angle
        'flip2';                %refocus pulse angle
		'p1';                   %excitation pulse duration
		'p2';                   %refocus pulse duration
        'pi';                   %inversion pi pulse
        'p1pat';                %excitation pulse shape
		'p2pat';                %refocus pulse shape
        'pipat';                %inversion pulse shape
        'epi_pc';               %'off', 'pointwise', 'linear', as described in feature list
        'epi_rev';              %1 reverses for odd echoes of epi, 0 does not reverse
        'rcvrs';                %number of y entries is number off channels
        'ne';                   %number of echoes
        'ns';                   %number of slices for standard experiment
        'ni';                   %number of images (?)
        'nseg';                 %number of segments or shots in multi-shot experiments
        'nblocks';              %number of blocks for arrayed experiments
        'nblk';                 %number of blocks in blockes, triggered acuisition (gems_x)
        'nsblock';              %number of slices per block (gems_x)
        'array';                %name(s) for the arrayed variables
        'arraydim';             %array dimension, should be nblocks from procpar but is not, though it is in the header of the fid
        'nt';                   %number of averages
        'np';                   %number of time domain points, divide by 2 because (real imaginary) pair counts as 2
        'nv';                   %number of phase encoded views = number of phase encodes (?)
        'fn';				    %number of readout points after zerofill (divide by four to get # complex points)
        'fn1';				    %number of slice encode points after zerofill
        'fn2';				    %number of slice encoded points after zerofill
        'fract_kx';				%number of acquired points past center in readout
        'fract_ky';             %number of acquired points past center in phase encode direction
		'dimX';                 %defines readout phase and slice
		'dimY';					%defines readout phase and slice
		'dimZ';					%defines readout phase and slice
        'ppe';                  %position phase direction
        'pro';                  %position readout direction
        'ppe2';                 %position second phase encode direction
        'pss';                  %slice positions for compressed experiment
        'pss0';                 %slice position offset
		'lro';                  %read out direction field of view
        'lpe';                  %phase direction field of view for shifting
        'lpe2';                 %phase direction field of view for shifting
        'thk';                  %slice thickness
        'nseg';                 %number of segments or shots in multishot experiment
        'etl';                  %echo train length
        'ms_intlv';             %interleaved: 1         sequential: 0
        'image';                %array of 0,1,-1 or -2 for epi; indicates reference scans
        'images';              %number of images if there are repeats
        'array';                %which parameters are arrayed
        'dcrmv';                %remove dc: 'y';        don't remove dc: 'n'    
        'nv2';                  %number of slice encoded views for 3D experiment
        'profile';              %'y','n' for 2D, 'yy', 'yn', 'ny', 'nn' for 3D experiment
        'altecho_reverse';      %'yy', 'yn', 'ny', 'nn' for readout/phase directions respectively
        'navigator';            %'y' if data was acquired with navigator echoes, 'n' if not
        'nav_echo';             %positions of navigator echoes in the echo train
        'navnp';                %number of points in the navigator echoes
        'nav_type';             %off, pointwise, linear, pairwise
        'pelist';               %Array of phasencodevalues
        'petable';              %external file containing array of phase encode values
        'ftproc';               %'y' apply weighting, 'n' do not apply
        'rcvrout';              %'i' output raw, pahse and image data for individual channels
        'raw';                  %'m','p',or'b' magnitude, phase or both raw
        'dmg';                  %'pa' output phase angle images
        'recon_force';          %1=reconstruct current data during scan
        'ppe2';                 %slice encoded direction shift
        'lpe2';                 %slice encoded field of view
        'comment';              %comment added by operator
        'sat';                  %satband check
        'seqfil';               %name of the sequencefile
        'pslabel';              %label of the pulse sequence
		'dro';					%diffusion multiplier dro
		'dpe';					%diffusion multiplier dpe
		'dsl';					%diffusion multiplier dsl
		'tdelta';				%spin echo little delta
		'tDELTA';               %Big DELTA for spin echo diffusion
        'TS';                   %storage time for the lig_stems sequence
        'ts';                   %storage time general
        'tabscheme';            %lig_stemems 2 lines per echo
        'epiflag';              %lig_stemems 2 lines per echo in epi mode
        'peshiftflag';          %lig_stemems shifted the k acquisition multiplier vpe
        'lambda';               %displacement encoding wave vector
		'max_bval';             %from the panel
        'fatoffs';              %fat offset from the panel
        'rflag';                %lag between slice selective and encoding gradient
        'acqlag';               %lag the acquisition to put kzero into the center of the at window
        'm1scale';              %scale the first moment
		'dixonf';               %dixon flag for fat separation, in UMS sequence
		'diffgrads';            %primary diffusion gradient flag in UMS se-sequence
		'diffgrads2';           %secondary diffusion gradient flag in UMS se-sequence
		'diffaxis';             %direction of primary diffusion gradient 
		'diffaxis2';            %direction of secondary diffusion gradient'
        'btarget';              %diffusion strength
        'Dtarget';              %diffusion strength
        'te_delay1';            %delay
        'te_delay2';            %delay
        'phi';                  %slice orientation
        'psi';                  %slice orientation
        'theta';                %slice orientation
        'svblist'};             %array of parameters to include in fdf headers
    
    %read the procparfile
    fid=fopen(procparfile);
    ii=0;
    while ~feof(fid);
        ii=ii+1;
        L{ii}=fgetl(fid);
    end
    fclose(fid);
    
    %consider the first line to identify the parameter line
    firstword=strtok(L);
    for jj=1:numel(parlist);
        line_number=find(strcmp(parlist{jj},firstword));
        if ~isempty(line_number);
            next_line=line_number+1;
            [dummy, remainder]=strtok(L{next_line});
            %if the remainder is a string, then it is enclosed in double
            %quotes; if it is numeric, no quotes
            quotes=regexp(remainder,'"');
            if isempty(quotes); %numeric
                out.(parlist{jj})=str2num(remainder);
            else
                out.(parlist{jj})=remainder(quotes(1)+1:quotes(2)-1);
            end
        end
    end
    
