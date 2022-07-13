function Uij = calculateUij(Fij)

maxTime = length(Fij);
Uij = cell(maxTime,1);

for i = 1:maxTime
    tic
  Uij{i} = funCalculateUij(Fij{i});
  disp(['Uij for time step ' num2str(i) ' complete']); 
    toc
end

end

function Uij = funCalculateUij(Fij)
Cij = cell(3,3);
Uij = cell(3,3);
Umat = cell(length(Fij{3,3}(:)),1);

Cij{1,1} = (Fij{1,1}.*Fij{1,1} + Fij{2,1}.*Fij{2,1} + Fij{3,1}.*Fij{3,1});
Cij{1,2} = (Fij{1,1}.*Fij{1,2} + Fij{2,1}.*Fij{2,2} + Fij{3,1}.*Fij{3,2});
Cij{1,3} = (Fij{1,1}.*Fij{1,3} + Fij{2,1}.*Fij{2,3} + Fij{3,1}.*Fij{3,3});
Cij{2,2} = (Fij{1,2}.*Fij{1,2} + Fij{2,2}.*Fij{2,2} + Fij{3,2}.*Fij{3,2});
Cij{2,3} = (Fij{1,2}.*Fij{1,3} + Fij{2,2}.*Fij{2,3} + Fij{3,2}.*Fij{3,3});
Cij{3,3} = (Fij{1,3}.*Fij{1,3} + Fij{2,3}.*Fij{2,3} + Fij{3,3}.*Fij{3,3});

for j=length(Fij{3,3}(:)):-1:1
    if ~isnan(Cij{1,1}(j))
        Umat{j} = sqrtm([Cij{1,1}(j), Cij{1,2}(j), Cij{1,3}(j); Cij{1,2}(j), Cij{2,2}(j), Cij{2,3}(j); Cij{1,3}(j), Cij{2,3}(j), Cij{3,3}(j)]);
    else
        Umat{j} = nan(3,3);
    end
        Uij{1,1}{j} = Umat{j}(1,1);
        Uij{1,2}{j} = Umat{j}(1,2);
        Uij{1,3}{j} = Umat{j}(1,3);
        Uij{2,2}{j} = Umat{j}(2,2);
        Uij{2,3}{j} = Umat{j}(2,3);
        Uij{3,3}{j} = Umat{j}(3,3);
end

Uij{1,1} = reshape(cell2mat(Uij{1,1}),size(Fij{1,1}));
Uij{1,2} = reshape(cell2mat(Uij{1,2}),size(Fij{1,1}));
Uij{1,3} = reshape(cell2mat(Uij{1,3}),size(Fij{1,1}));
Uij{2,2} = reshape(cell2mat(Uij{2,2}),size(Fij{1,1}));
Uij{2,3} = reshape(cell2mat(Uij{2,3}),size(Fij{1,1}));
Uij{3,3} = reshape(cell2mat(Uij{3,3}),size(Fij{1,1}));
%Umat = sqrtm([Cij{1,1}, Cij{1,2}, Cij{1,3}; Cij{1,2}, Cij{2,2}, Cij{2,3}; Cij{1,3}, Cij{2,3}, Cij{3,3}]);

end


