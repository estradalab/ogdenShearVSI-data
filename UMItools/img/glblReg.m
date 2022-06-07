function glbl = glblReg(rootname)

imgs = read_img_series(rootname);
glbl = mean(imgs,2);
save glbalMean.dat glbl -ASCII

plot(glbl)

return
