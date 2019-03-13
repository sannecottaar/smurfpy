
#
# Script to Count number of files in 'auto_removed' folders

import glob,os
import numpy as np
import matplotlib.pyplot as plt
from obspy import read
import shutil

# Find list of stations directories
stations = glob.glob('Data/*')
print ('number of stations', len(stations))

# Directories to remove
remove_dirs = ['Processed_originals','Travel_time_added','auto_removed']
#
# Script to Count number of files not in 'Originals' folders

# Find list of stations directories
stations = glob.glob('Data/*')

# Loop through stations

count_removed = 0

for stadir in stations:
        print(stadir)
        for dr in remove_dirs:
                stalist = glob.glob(stadir + '/'+dr+'/*PICKLE')
                for st in stalist:
                        os.remove(st)
                        count_removed+=1


print ("Number discarded: ", count_removed)
