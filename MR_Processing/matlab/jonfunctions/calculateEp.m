function [Q,En] = calculateEp(Eij)
maxTime = length(Eij);
En = cell(maxTime,1);
Q = cell(maxTime,1);

for i = 1:maxTime
    tic
  [Q{i},En{i}] = funCalculateEp(Eij{i});  
  disp(['En for time step ' num2str(i) ' complete']); 
    toc
end

end

function [Q,En] = funCalculateEp(Eij)

Q = cell(3,3);
En = cell(3,3);
Emat = cell(length(Eij{3,3}(:)),1);
%Qmat = cell(length(Eij{3,3}(:)),1);

for j=1:length(Eij{3,3}(:));
    
Emat{j} = [Eij{1,1}(j), Eij{1,2}(j), Eij{1,3}(j); 
    Eij{1,2}(j), Eij{2,2}(j), Eij{2,3}(j); 
    Eij{1,3}(j), Eij{2,3}(j), Eij{3,3}(j)];

[Qmat, Enmat] = eig(Emat{j});

Q{1,1}{j} = Qmat(1,1);
Q{1,2}{j} = Qmat(1,2);
Q{1,3}{j} = Qmat(1,3);
Q{2,1}{j} = Qmat(2,1);
Q{2,2}{j} = Qmat(2,2);
Q{2,3}{j} = Qmat(2,3);
Q{3,1}{j} = Qmat(3,1);
Q{3,2}{j} = Qmat(3,2);
Q{3,3}{j} = Qmat(3,3);

En{1}{j} = Enmat(1,1); 
En{2}{j} = Enmat(2,2);
En{3}{j} = Enmat(3,3);
end

Q{1,1} = reshape(cell2mat(Q{1,1}),size(Eij{1,1}));
Q{1,2} = reshape(cell2mat(Q{1,2}),size(Eij{1,1}));
Q{1,3} = reshape(cell2mat(Q{1,3}),size(Eij{1,1}));
Q{2,1} = reshape(cell2mat(Q{2,2}),size(Eij{1,1}));
Q{2,2} = reshape(cell2mat(Q{2,2}),size(Eij{1,1}));
Q{2,3} = reshape(cell2mat(Q{2,3}),size(Eij{1,1}));
Q{3,1} = reshape(cell2mat(Q{3,1}),size(Eij{1,1}));
Q{3,2} = reshape(cell2mat(Q{3,2}),size(Eij{1,1}));
Q{3,3} = reshape(cell2mat(Q{3,3}),size(Eij{1,1}));
En{1} = reshape(cell2mat(En{1}),size(Eij{1,1}));
En{2} = reshape(cell2mat(En{2}),size(Eij{1,1}));
En{3} = reshape(cell2mat(En{3}),size(Eij{1,1}));

end