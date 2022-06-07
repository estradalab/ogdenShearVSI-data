function result = rightleft3(data, plane)
% function result = rightleft3(data, plane)
%
% Luis Hernandez
% Last Edit: 7-3-98
%
% This is like rightleft - No plotting anything, though.
% given a set of 3D points and the 3D coordinates of a plane
% returns how many are on the right (+x) and how many are on the
% left of the plane
% 
% The function also computes the weighted sum of the voxel locations
% of the two sides of the plane and returns a matrix of the form
%
% [x_right y_right z_right total_right_weight rightcount;
%  x_left  y_left  z_left total_left_weight leftcount]
%

%disp('starting rightleft3');

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
 
 %hold off
 
 % Compute the weighted sums
 ws_right = w_sum(rightdata) ;
 ws_left =  w_sum(leftdata);
 
 szr = size(rightdata);
 szl = size(leftdata);
 
 result = [ws_right szr(1);
           ws_left szl(1)];
 
% disp('ending rightleft3');
 
 return
 
 %%%%%%%%%%%%%%%%%%%%%%%%


function result = w_sum(data)
%
% computes the weighted average position from the points 
% weighted by their intensity
%

weighted = zeros(1,4);

sz = size(data);

   for i=1:sz(1)
      weighted(1:3) = weighted(1:3) + data(i,1:3) * data(i,4);
      weighted(1:4) = weighted(1:4) + data(i,1:4);
   end
   
	
   result = weighted;
   
 return 
