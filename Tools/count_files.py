
#
# Script to Count number of files in 'auto_removed' folders

import glob
import numpy as np

# Find list of stations directories
stations = glob.glob('../Data/*')
print ('number of stations', len(stations))

station_names = []
station_data = []

# Script to Count number of files not in 'Originals' folders

# Loop through stations

count_total = 0
count_removed = 0

for stadir in stations:
        stalist = glob.glob(stadir + '/*PICKLE')
        count_total = count_total + len(stalist)
        stalist = glob.glob(stadir + '/auto_removed/*PICKLE')
        count_removed = count_removed + len(stalist)

print ("Number of seismograms: ", count_total)
print ("Number discarded: ", count_removed)
