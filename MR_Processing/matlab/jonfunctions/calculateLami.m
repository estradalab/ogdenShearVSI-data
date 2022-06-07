function Lami = calculateLami(Uij)

maxTime = length(Uij);
Lami = cell(maxTime,1);

for i = 1:maxTime
    tic
  Lami{i} = funCalculateLami(Uij{i});
  disp(['Lambda_i for time step ' num2str(i) ' complete']); 
    toc
end

end


function Lami = funCalculateLami(Uij)
Lami = cell(3,1);

for j=length(Uij{3,3}(:)):-1:1
    if ~isnan(Uij{1,1}(j))
        [~,Lmat{j}] = eigenshuffle([Uij{1,1}(j), Uij{1,2}(j), Uij{1,3}(j); Uij{1,2}(j), Uij{2,2}(j), Uij{2,3}(j); Uij{1,3}(j), Uij{2,3}(j), Uij{3,3}(j)]);
    else
        Lmat{j} = nan(3,3);
    end
        Lami{1}{j} = Lmat{j}(1);
        Lami{2}{j} = Lmat{j}(2);
        Lami{3}{j} = Lmat{j}(3);
end

Lami{1} = reshape(cell2mat(Lami{1}),size(Uij{1,1}));
Lami{2} = reshape(cell2mat(Lami{2}),size(Uij{1,1}));
Lami{3} = reshape(cell2mat(Lami{3}),size(Uij{1,1}));

end