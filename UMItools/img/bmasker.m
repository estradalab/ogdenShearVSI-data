function bmasker(root, thres)
% function bmasker(root, thres)
% 
% makes a binary mask at a specified lower threshold.
% and writes it to the file : bmask.img
%

[img h] = read_img2(root);

bmask = img;
bmask(img<thres) = 0;
bmask(img>=thres) = 1;

write_img('bmask.img', bmask, h);

lightbox(bmask);
return


