function [freqaxis,ft,maxfrequency]=varianspuls(filename)
%[freqaxis,ft,maxfrequency]=varianspuls(filename)
%open the fid, fft it, get the peak frequency;

ft=[];
frequency=[];

if isdir(filename);
    %1. check if the directory is indeed an spuls directory
    if ~isempty(strfind(filename,'spuls'));
        opar=getparameters([filename '/procpar']);
        filename=[filename '/fid'];
    else
        display(['Error using varianspuls: ' filename ' is not an spuls directory.']);
        return
    end
else
    %2. a file name was handed down
    [PATH,NAME,EXT]=fileparts(filename);
    is_spuls=false; %default assumption
    if strcmp(NAME,'fid'); %it is an fid
        if isempty(PATH);
            opar=getparameters('procpar');
          
        else
            opar=getparameters([PATH '/procpar']);
        end
        is_spuls=strcmp(opar.seqfil,'spuls');
    end
    if ~is_spuls;
        return
    end
end
    

file_id=fopen(filename, 'r', 'b');          %always big endian format, even for Linux

%Start reading the file header to get the first 6 elements
[hdr, count] = fread(file_id, 6, 'int32'); %count is number of

header.nblocks = hdr(1);        %data blocks in a file (second dimension, if arrayed acquisition)
header.ntraces = hdr(2);        %number of traces
header.complex_td= hdr(3)/2;    %real and imaginary points (first dimension = time);
ebytes = hdr(4);                % 4 bytes per point	(double precision real, double precision imag)
tbytes = hdr(5);                %np*ebytes
bbytes = hdr(6);                %ntraces*tbytes+nhead*sizeof(datablockhead=28bytes)-total

%Read 2 16bit integers to see the data type and the status of the binary
[data_holder, count] = fread(file_id, 2, 'int16');
vers_file_id = data_holder(1);
status = data_holder(2);            %status determines whether data are store as 16 or 32 bit integers or floats
BinStat = fliplr(dec2bin(status));  %determine the binary format from the status info
%xx00xxxx means 16 bit data
%xx10xxxx means 32 bit data
%else means float

[data_holder, count] = fread(file_id, 1, 'int32');
nhead = data_holder(1);

data_holder = zeros(hdr(3)*hdr(2), 1);

%INT16 data type **********************************************************
if ((BinStat(3) == '0') & (BinStat(4) == '0'))
    for nblocks_count=1:header.nblocks
        for nhead_count=1:nhead
            fread(file_id, 14, 'int16');			%read in the block headers (28 bytes)
        end
        [b, count] = fread(file_id, hdr(2)*hdr(3), 'int16');   	%read in the actual data (ntraces*np)
        %put the data in the vector
        data_holder( ((nblocks_count-1)*hdr(2)*hdr(3)+1) : nblocks_count*hdr(2)*hdr(3) ) = b;
    end
    %INT32 data type **********************************************************
elseif ((BinStat(3) == '1') & (BinStat(4) == '0'))
    for nblocks_count=1:header.nblocks
        for nhead_count=1:nhead
            [d,count]=fread(file_id, 14, 'int16');			%read in the block headers (28 bytes)
            bhead=d;
        end
        [b, count] = fread(file_id, hdr(2)*hdr(3), 'int32');   	%read in the actual data (ntraces*np)
        %put the data in the vector
        data_holder( ((nblocks_count-1)*hdr(2)*hdr(3)+1) : nblocks_counthdr(2)*hdr(3) ) = b;
    end
    %FLOAT data type **********************************************************
else
    for nblocks_count=1:header.nblocks
        for nhead_count=1:nhead
            fread(file_id, 14, 'int16');			%read in the block headers (28 bytes)
        end
        [b, count] = fread(file_id, hdr(2)*hdr(3), 'float');   	%read in the actual data (ntraces*np)
        %put the data in the vector
        data_holder( ((nblocks_count-1)*hdr(2)*hdr(3)+1) : nblocks_count*hdr(2)*hdr(3) ) = b;
    end
end

%real imaginary pairs, create complex array
acquireddata=complex(data_holder(1:2:end),data_holder(2:2:end));
acquireddata=reshape(acquireddata,[header.complex_td numel(acquireddata)/header.complex_td]);

ft=fftshift(fft(acquireddata));
nt=numel(acquireddata);
dwell=opar.at/nt;
freqrange=1/dwell;
freqaxis=[nt/2:-1:-nt/2+1]/nt*freqrange;
[dummy,maxind]=max(abs(ft));
maxfrequency=freqaxis(maxind);

