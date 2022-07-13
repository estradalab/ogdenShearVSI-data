function out=matrix_satrec_polyfit(varargin)
%out=matrix_satrec_polyfit(indata,taxis[,options?])
%fits the recovery vs time to linear plus quadratic terms
%to get the ratio (spindensity/Tau) for each tissue.


indata=varargin{1};
taxis=varargin{2};

if ndims(indata) == 3;                          %one slice, multiple times
        out = matrix_satrec_polyfit_3d(indata,taxis);
elseif ndims(indata) ==4;                       %multislice, multiple times
    [ntd nph nsl nti] = size(indata);
    out.amplitude=zeros(ntd,nph,nsl);
    out.Tau=zeros(ntd,nph,nsl);
    
    for ns=1:nsl;
        oout= matrix_satrec_polyfit_3d(squeeze(indata(:,:,ns,:)),taxis);
        out.A(:,:,ns)=oout.A(:,:);
        out.B(:,:,ns)=oout.B(:,:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saturation recovery poly fit linear plus quadratic term:
function out=matrix_satrec_polyfit_3d(varargin)


incube=varargin{1};
taxis=varargin{2};
[tmax,tmaxind]=max(taxis);
[tmin,tminind]=min(taxis);

sincube=size(incube);

if ndims(incube) == 3;
    %third dimension is the echotime dimension
    phasemat=incube(:,:,tmaxind)./abs(incube(:,:,tmaxind)); %phase of the image with the shortest wait time
    phincube=zeros(size(incube));
    
    %1. rotate the phase
    for nt=1:numel(taxis);
         phincube(:,:,nt)= +incube(:,:,nt)./phasemat;
    end
    
    %2.calculate the coefficients of the time matrix
    for jj=1:4;
        tpow(jj)=sum(taxis.^jj);
    end
    invT=1/(tpow(2)*tpow(4)-tpow(3)^2)*[tpow(4) -tpow(3); -tpow(3) tpow(2)];
    
    %3. calculate the signal time power of t coefficients
    S1=zeros(sincube(1:2));
    S2=zeros(sincube(1:2));
    for jj=1:numel(taxis);
        S1=S1+abs(incube(:,:,jj))*taxis(jj);
        S2=S2+abs(incube(:,:,jj))*taxis(jj)^2;
    end
        
    %4. calculate the coefficients
    A=invT(1,1)*S1+invT(1,2)*S2;
    B=invT(2,1)*S1+invT(2,2)*S2;
    
end

out.A=A;
out.B=B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

