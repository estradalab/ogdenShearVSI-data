% Data
x=[0:5:40];
y= 3*x + x.^3 - 4*x.^5 + sin(x/5) + 1.2.*cos(x/3);
xsize = size(x,2);
speriod = x(2)-x(1);


plot(x,y);
hold on


% xcoordinate to be interpolated
for xint= 0:1:40;
   wtsize = 3;
   xwt=zeros(1,wtsize*2);
   ywt=xwt;


   if ~isempty( find(x==xint) )
      yint = y(find( x == xint));
   else
   
      lindices = find( x > xint - wtsize*speriod   &   x < xint)
      rindices = find( x < xint + wtsize*speriod   &   x > xint)
      shift = size(rindices,2) - size(lindices,2)
      
      indices = find( x < xint + (wtsize + shift)*speriod & ...
            					x > xint - (wtsize - shift)*speriod )
            
      for i= 1:2*wtsize
         xwt(i) =  x(indices(i) ) - xint  ;
      end
   
      ywt = sin( (pi/speriod)* xwt ) ./ ( (pi/speriod)* xwt  )
      yint = sum(ywt .* y(indices(:))  ) / sum(ywt)
   
   end
   
   plot(xint,yint,'g*');
   
end
