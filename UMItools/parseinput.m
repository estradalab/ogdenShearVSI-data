function parseinput(inputwordlist,varargin_inputlist )
%parseinput(inputputwords,cellofinputs )
%   inputwords is a cell array of strings, to be matched with the cell varargin_inputslist

num_words  =numel(inputwordlist);
num_inputs =numel(varargin_inputlist);

%set the defaults
for nw=1:num_words;
     assignin('caller',[inputwordlist{nw} '_flag'],false);
end


for nw=1:num_words;
    for ni=1:num_inputs;
        if strcmp(inputwordlist{nw},varargin_inputlist{ni});
            assignin('caller',[inputwordlist{nw} '_flag'],true);
            
            if ni<num_inputs && isnumeric(varargin_inputlist{ni+1});
            assignin('caller',[inputwordlist{nw} '_value'],varargin_inputlist{ni+1});
            end
        end
    end
end

