#!/bin/bash

# 1. copy the nmr data of the mirror on enthusiasm 
rsync --update -avz uscheven@enthusiasm.engin.umich.edu:/mirror/uscheven/vnmrsys/data/Ligaments Documents/UMI/vnmrsys/data/Ligaments

# 1b. pull the desktop
# rsync -avz uscheven@enthusiasm.engin.umich.edu:/mirror/uscheven/Desktop/ Documents/UMI/Desktop_varian

# 2. pull the xugroup data folder
# rsync -avz uscheven@enthusiasm.engin.umich.edu:/mirror/xugroup/vnmrsys/data/ Documents/UMI/vnmrsys/data/xugroup
