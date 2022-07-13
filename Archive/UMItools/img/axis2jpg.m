function axis2jpg(name)

f = getframe;
[i map] = frame2im(f);

imwrite(i,map,name,'jpg');