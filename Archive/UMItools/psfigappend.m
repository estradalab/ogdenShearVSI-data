function out=psfigappend(handlearray, filename)
% out=psfigappend(currentfigurehandle, filename);
% appends the figure specified by the input handle to the pdf file
% specified in name. if the file does not exist, cretes it


numfigs=numel(handlearray);

for jj=1:numfigs;
    figurehandle=mandlearray(jj);
    set(figurehandle,'PaperOrientation','landscape','paperpositionmode','manual','paperposition', [0.25 0.25 10.5 7.5]);
    
    %print the figure to a tempfile
    tempfilename=['temp_pdf' num2str(jj) '.pdf'];
    
    print(figurehandle, '-dpdf', tempfilename);
end

%concatenate the pdf slides