function [Rij, Uij, angles] = calculateRij(Fij)

maxTime = length(Fij);
Rij = cell(maxTime,1);
Uij = cell(maxTime,1);
angles = cell(maxTime,1);

for i = 1:maxTime
    tic
  Uij{i} = funCalculateUij(Fij{i});
  [Rij{i},angles{i}] = funCalculateRij(Fij{i}, Uij{i});  
  disp(['Rij for time step ' num2str(i) ' complete']); 
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

function [Rij,angles] = funCalculateRij(Fij, Uij)
Rij = cell(3,3);
Umat = cell(length(Fij{3,3}(:)),1);
Rmat = cell(length(Fij{3,3}(:)),1);
Fmat = cell(length(Fij{3,3}(:)),1);
angles = cell(3,1);


for j=1:length(Fij{3,3}(:));
Umat{j} = [Uij{1,1}(j), Uij{1,2}(j), Uij{1,3}(j); Uij{1,2}(j), Uij{2,2}(j), Uij{2,3}(j); Uij{1,3}(j), Uij{2,3}(j), Uij{3,3}(j)];
Fmat{j} = [Fij{1,1}(j), Fij{1,2}(j), Fij{1,3}(j); Fij{2,1}(j), Fij{2,2}(j), Fij{2,3}(j); Fij{3,1}(j), Fij{3,2}(j), Fij{3,3}(j)];
Rmat{j} = Fmat{j}/Umat{j};
Rij{1,1}{j} = Rmat{j}(1,1);
Rij{1,2}{j} = Rmat{j}(1,2);
Rij{1,3}{j} = Rmat{j}(1,3);
Rij{2,1}{j} = Rmat{j}(2,1);
Rij{2,2}{j} = Rmat{j}(2,2);
Rij{2,3}{j} = Rmat{j}(2,3);
Rij{3,1}{j} = Rmat{j}(3,1);
Rij{3,2}{j} = Rmat{j}(3,2);
Rij{3,3}{j} = Rmat{j}(3,3);

angles{1}{j}=atan2(Rij{3,2}{j},Rij{3,3}{j})*180/pi; %range: -pi to pi
angles{2}{j}=atan2(-Rij{3,1}{j},sqrt(Rij{3,2}{j}.^2+Rij{3,3}{j}.^2))*180/pi; %range: -pi/2 to pi/2
angles{3}{j}=atan2(Rij{2,1}{j},Rij{1,1}{j})*180/pi; %range: -pi to pi
end

Rij{1,1} = reshape(cell2mat(Rij{1,1}),size(Fij{1,1}));
Rij{1,2} = reshape(cell2mat(Rij{1,2}),size(Fij{1,1}));
Rij{1,3} = reshape(cell2mat(Rij{1,3}),size(Fij{1,1}));
Rij{2,1} = reshape(cell2mat(Rij{2,2}),size(Fij{1,1}));
Rij{2,2} = reshape(cell2mat(Rij{2,2}),size(Fij{1,1}));
Rij{2,3} = reshape(cell2mat(Rij{2,3}),size(Fij{1,1}));
Rij{3,1} = reshape(cell2mat(Rij{3,1}),size(Fij{1,1}));
Rij{3,2} = reshape(cell2mat(Rij{3,2}),size(Fij{1,1}));
Rij{3,3} = reshape(cell2mat(Rij{3,3}),size(Fij{1,1}));
angles{1} = reshape(cell2mat(angles{1}),size(Fij{1,1}));
angles{2} = reshape(cell2mat(angles{2}),size(Fij{1,1}));
angles{3} = reshape(cell2mat(angles{3}),size(Fij{1,1}));
end

