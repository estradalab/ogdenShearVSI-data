%function build_img name

	
	fmt = 'uint16';
	[pFile,messg] = fopen(name, 'r');
	if pFile == -1
		disp(messg);   
		return;
	end
  

	dummy = fread(pFile, 80,fmt);
	data = fread(pFile,[128 128], fmt) ;	
	fclose(pFile);
	
	colormap(gray);
	imagesc(data);


