function set_func_colors2

    mygray = [0:1/255:1]' * [1 1 1]; 
    
	myhot = [0:3:255]' * [1 0 0] ;
	tmp =   [0:3:255]' * [0 1 0] ;
	myhot = [myhot; tmp];
	tmp =   [0:3:255]' * [0 0 1];
	myhot =  [myhot;  tmp]/256;

	myhot(round(256/3): 256, 1) = 1;
	myhot(round(256*2/3):256,2) = 1;

	myblue = myhot;
	myblue(:,1) = myhot(:,3);
	myblue(:,3) = myhot(:,1);

	mygreen = myhot;
	mygreen(:,3) = 0;%mygray(:,1);
	mygreen(:,1) = 0;
	mygreen(:,2) = 1;

 	mymap = [mygray; myhot ; myblue ; mygreen];
    colormap(mymap)