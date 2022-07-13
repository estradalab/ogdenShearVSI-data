function h=plotms(out,varargin)
%h=plotms(instructure,varargin)
%function to plot the images produced by varianms, or a matrix
%instructure is the structure returned by varianms
%options:
%'clim', [cmin cmax],   %sets the color scale
%'pix',             %x and y axes in pixels
%'layout', [nrows ncols]
%'square' %square layout;

options=varargin;
if isstruct(out);
    %make axes for plotting
    si=size(out.image);
    if isfield(out,'pars');
        if out.pars.phi==0;
            %RO is ydirection, PE is x direction
            out.pars.yaxis=-((1:si(1))/si(1)*out.pars.lro+out.pars.pro)*10;
            out.pars.xaxis=((1:si(2))/si(2)*out.pars.lpe+out.pars.ppe)*10;
        elseif out.pars.phi==90;
            out.pars.yaxis=-((1:si(1))/si(1)*out.pars.lpe+out.pars.ppe)*10;
            out.pars.xaxis=((1:si(2))/si(2)*out.pars.lro+out.pars.pro)*10;
        else
            out.pars.yaxis=-((1:si(1))/si(1)*out.pars.lro)*10;
            out.pars.xaxis=((1:si(2))/si(2)*out.pars.lpe)*10;
        end
        h=slice_surveyplot(out, options{:});
    end
else
    %the input variable out is a matrix or nD double
    dummy.image=out;
    h=slice_surveyplot(dummy, options{:});
end

