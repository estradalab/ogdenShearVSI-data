function pout=p_semgems(varargin)
%plots amplitude on phase for all echoes, one figure per slice;

in=varargin{1};
if nargin>1 && isnumeric(varargin{2});
	sli=varargin{2};
	nsl=numel(sli);
else
	if isfield(in.pars,'ns'); nsl=in.pars.ns; end	
	if isfield(in.pars,'slices'); nsl=in.pars.slices; end
	sli=1:nsl;
end

if isfield(in.pars,'ne'); nec=in.pars.ne; end
if isfield(in.pars,'echoes'); nec=in.pars.ne; end
if isfield(in.pars,'plotlim'); 
    plimflag=true;
    clim = in.pars.plotlim;
else
    plimflag=false;
end


if any(strcmp(varargin,'bw'));
    cm=gray(64);
else
    cm=jet(64);
end

%optional plot;
npx=250; % pixels per patch
maxamp=abs(max(in.image(:)));
if nec>1 && nsl>1; %multi-slice multi-echo
	for ns=1:nsl;
		figure('Name', ['slice ' num2str(sli(ns))],'Position',[50 50 nec*npx, 2*npx],'colormap',cm);

		for ne=1:nec;
			%magnitude
			subplot(2,nec,ne);
			imagesc(abs(in.image(:,:,sli(ns),ne)));
            if ~plimflag;
                set(gca,'clim',[0 maxamp]);
            else
                set(gca,'clim',clim);
            end
            
			title(['abs, echo ' num2str(ne)]);
			axis square;
			%phase
			subplot(2,nec,ne+nec);
			imagesc(angle(in.image(:,:,sli(ns),ne)));
			set(gca,'clim',[-pi pi]);
			axis square;
			title(['phase, echo ' num2str(ne)]);
		end
	end
end

%single echo, multislice
if nec==1 && nsl>1;
    for ns=1:nsl;
        figure('Name', ['slice ' num2str(sli(ns))],'Position',[50 50 2*npx, npx]);
        %magnitude
        subplot(1,2,1);
        imagesc(abs(in.image(:,:,sli(ns))));
        if ~plimflag;
            set(gca,'clim',[0 maxamp]);
        else
            set(gca,'clim',clim);
        end
        title('abs');
        axis square;
        %phase
        subplot(1,2,2);
        imagesc(angle(in.image(:,:,sli(ns))));
        set(gca,'clim',[-pi pi]);
        axis square;
        title('phase');
    end
end

%single slice multi-echo
if nec>1 && nsl==1;
    ns=1;
    figure('Name', ['slice ' num2str(sli(ns))],'Position',[50 50 nec*npx, 2*npx]);
    for ne=1:nec;
        %magnitude
        subplot(2,nec,ne);
        if ndims(in.image)==3;
            imagesc(abs(in.image(:,:,ne)));
        else
            imagesc(abs(in.image(:,:,sli(ns),ne)));
        end
        set(gca,'clim',[0 maxamp]);
        title(['abs, echo ' num2str(ne)]);
        axis square;
        %phase
        subplot(2,nec,ne+nec);
        if ndims(in.image)==3;
            imagesc(angle(in.image(:,:,ne)));
        else
            imagesc(angle(in.image(:,:,sli(ns),ne)));
        end
        set(gca,'clim',[-pi pi]);
        axis square;
        title(['phase, echo ' num2str(ne)]);
    end
end
