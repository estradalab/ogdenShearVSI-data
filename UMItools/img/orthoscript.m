  args.ROIsize = 0;      
      args.ROItype = 'maskFile'; 
      args.threshold = 1;
      args.onsets = [8:8:152];
      args.window = 10;
      args.spm_file = '/raid/data6/ROIs/resliced_masks/rLR_BA17_d3.img';
      args.spm_file2 = [];
      args.anat_file = '/usr/local/spm99/templates/T1.img';
      args.tseries_path = '/raid/data6/EQdiff/030910er/run_1/'
      args.tseries_file =  '/raid/data6/EQdiff/030910er/run_1/snravol_e2231_09_10_103_0160.img';
      args.tseries_file2 = [];
      args.doDetrend = 0;
      args.doGfilter = 0;
      args.doFFT = 1;
      args.ignore_origin = 0;
      args.wscale = [];
      args.interact = 0;
      args.xyz=[];
      args.mask_file = '/raid/data6/ROIs/resliced_masks/rLR_BA17_d3.img';
      args.output_name = 'Ortho';
      args.voxFile = [];

ortho2005(args)
