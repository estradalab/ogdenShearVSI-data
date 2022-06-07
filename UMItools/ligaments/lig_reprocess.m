function mssg = lig_reprocess
% looks for matlab workspaces  DESTE_strains_####_#####_#####_#####.mat
%extracts the time stamps and re-runs lig_deste_3d on those

d=dir;

for jj=1:numel(d)
	
	if contains(d(jj).name,'DESTE_strains') && ~contains(d(jj).name,'.pdf')
		current_DESTE=(d(jj).name);
		
		underscores=regexp(current_DESTE,'_');
		
		for kk=2:numel(underscores); %the first underscore is between DESTE and strains
			timestamp{kk-1}=str2num(current_DESTE(underscores(kk)+(1:4)));
		end
		
		lig_deste_3d(timestamp{:});
	end
end

