function build_img(name)
% build_img(name)

	
	fmt = 'uint8';
	[pFile,messg] = fopen(name, 'r');
	if pFile == -1
		disp(messg);   
		return;
	end
  

	dummy = fread(pFile, 80,fmt);
	data = fread(pFile,[64 64], fmt) ;	


	colormap(gray);
	imagesc(data);

return
