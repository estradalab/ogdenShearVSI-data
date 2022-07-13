% sample script to run through a bunch of runs and subjects and average
% them together:

% This is where we set all the defaults
% change things to reflect your desired analysis
% this are the most important ones at the top.
args.onsets = [];
args.window = 10;
args.spm_file = []; %'/apollohome3/reward/nalt/spm210/pairedT/n7/spmT_0003.img';
args.anat_file = [];%'/matrixhome4/anw3049819mr.img';
args.mask_file = '';
args.ROItype = 'voxFile'; 
args.threshold = 3.5;
args.ROIsize = 0;
args.spm_file2 = [];%'/apollohome3/reward/nalt/spm210/pairedT/n7/spmT_0002.img';
args.tseries_path = '';
args.tseries_file = [];
args.tseries_file2 =[];
args.doDetrend = 0;
args.doGfilter = 0;
args.doFFT = 0;
args.ignore_origin = 0;
args.wscale = [];
args.interact = 0;
args.xyz=[];
args.output_name = 'Ortho';
args.voxFile = 'Ortho_voxels.dat';
args.doMovie = 0;
args.causalMap = 1;


addpath('/apollohome3/luistools/ortho2005')
addpath('/apollohome3/luistools/')
subj_dirs = {'050629zb', '050630mh'}
rootdir = '/apollohome3/reward/nalt'
runs = {'run_01'};% 'run_2' 'run_3' 'run_4' 'run_5' 'run_6' 'run_7' 'run_8' 'run_9'}

cd (rootdir)


for s = 1:length(subj_dirs)
    % get into every directory, one at a time
    
    cd (cell2mat(subj_dirs(s)))
    fprintf('\nDoing subject %s', cell2mat(subj_dirs(s)));
    avgtimeseries = zeros(args.window,1);
    DMfile = sprintf('%s/spmhrf/%s/pre/SPM.mat',rootdir,cell2mat(subj_dirs(s)))
    load(DMfile);
    dm = SPM.xX.X;
    der = [diff(dm(:,1)) ; 0];
    dm = [dm der];
    
    for r=1:length(runs)

        runstr = sprintf('func/%s/a_img', cell2mat(runs(r)))
        cd(runstr)
    	save DesMat.dat dm -ASCII
        
        fprintf('\nDoing run .... %s', cell2mat(runs(r)));
    	tseries = dir('ravol_*.img')
    	args.tseries_file = tseries(1).name;
        %args
    	ortho2005(args)
        %load Ortho_avg.dat
        %avgtimeseries = avgtimeseries + Ortho_avg(:,1);
        
        cd ../../..
        drawnow
    end
    cd ..
    avgtimeseries = avgtimeseries/length(runs);
    save RUN_avg.dat avgtimeseries -ASCII
    cd(rootdir)
    close all
end


fprintf('\All Done')
