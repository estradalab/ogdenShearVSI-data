function result = rightleft2( data, plane)
% function result = rightleft2( data, plane)
%
% Luis Hernandez
% Last Edit: 3-24-98
%
% This is like rightleft - No plotting anything, though.
% given a set of 3D points and the 3D coordinates of a plane
% returns how many are on the right (+x) and how many are on the
% left of the plane
% 
% The function also computes the centers of gravity 
% of the two sides of the plane by computing the weighted average
% position of the voxels.

disp('starting rightleft2');

global slice_number

% Determine two vectors in the plane.  Their cross-product is the
% normal of the plane

vector1 = plane(1,:) - plane(2,:);
vector2 = plane(3,:) - plane(2,:);
normal = cross(vector1, vector2);

% (xo, yo, zo)is the point that the two vectors have in common
xo = plane(2,1);
yo = plane(2,2);
zo = plane(2,3);


right = 0;
left = 0;

rightdata = [0 0 0 0];
leftdata = [0 0 0 0];

sz = size(data);
for i=1:sz(1)

   x = data(i,1);
   y = data(i,2);
   z = data(i,3);

   %hold on

   % equation of the plane is determined from a line in the plane
   % and the normal of the plane.
 
  if x > xo + ( normal(2)*(yo-y) + normal(3)*(zo-z) )/normal(1) 
     
     left = left+1;
     leftdata = [leftdata; data(i,:)];

     
   else
      right = right+1;    
      rightdata = [rightdata; data(i,:)];

   end

 end
% remove the top row (contains zeros)
sr = size(rightdata);
sl = size(leftdata);

rightdata = rightdata(2:sr(1), :);
leftdata = leftdata(2:sl(1), :);

% Plot the SPM data points
hold on
plot3(rightdata(:,1),rightdata(:,2),rightdata(:,3),'r*')
plot3(leftdata(:,1),leftdata(:,2),leftdata(:,3),'b*')
hold off

disp('-------- Right voxels, Left voxels------')
[right left]
result = [right left];

disp('--------Center of Mass---------------');
wc = w_ave(data)

disp('----Right, Left Center of Mass--------');
if ~isempty(rightdata)
   rightcenter = w_ave(rightdata)
end	

disp('----Right, Left Center of Mass--------')
if ~isempty(leftdata)
    leftcenter = w_ave(leftdata)
end

[fn,pn] = uiputfile('*.txt','Save Stats As ...');
fn = strcat(pn,fn);

fp = fopen(fn,'w');

fprintf(fp,'Right voxels:	\t%d\n',right);
fprintf(fp,'Left  voxels:	\t%d\n',left);
fprintf(fp,'Center of Mass: \t%f\t%f\t%f\n',wc);
fprintf(fp,'Right Centroid: \t%f\t%f\t%f\n',rightcenter);
fprintf(fp,'Left  Centroid: \t%f\t%f\t%f\n',leftcenter);

fclose (fp);
disp('ending rightleft2');
return

%%%%%%%%%%%%%%%%%%%%%%%%


function result = w_ave(data)

% computes the weighted average position from the points weighted by their intensity

	sz = size(data);
	for i=1:3
	   weighted(:,i) = data(:,i) .* data(:,4);
	end

	weighted_sum = [sum(weighted(:,1))  sum(weighted(:,2)) sum(weighted(:,3))];
	total_mass = sum(data(:,4));
	wa = weighted_sum / total_mass;

	%plot3(wa(:,1), wa(:,2), wa(:,3) , 'g*')

	result = wa;
return 
