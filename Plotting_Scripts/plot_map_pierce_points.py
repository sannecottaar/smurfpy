###########################
# Plot map with piercepoints
############################

import sys
import mpl_toolkits
import mpl_toolkits.basemap
from mpl_toolkits.basemap import Basemap

import numpy as np
# from matplotlib.mlab import griddata
import scipy
from scipy import interpolate
import matplotlib.pyplot as plt
import subprocess
import glob

# Set file with pierce points
depth = '410'   # depth of pierce points to plot
phase = 'P410s' # phase of pierce points to plot 
rffilter = 'jgf1'
piercelist = ['../Migration_Scripts/PP_'+depth+'km_'+phase+'_'+rffilter+'.txt']

# Read in pierce points
lonpp = []
latpp = []
print(piercelist)
for filename in piercelist:
    print(filename)
    rd = open(filename, 'r')

    for line in rd.readlines():
        val = line.split()
        lonpp.append(float(val[4]))
        latpp.append(float(val[3]))

    rd.close()


plt.figure(figsize=(6, 10))

# Set bounds of map based on piercepoints
lonmin = np.min(lonpp) - 2
lonmax = np.max(lonpp) + 2
latmin = np.min(latpp) - 2
latmax = np.max(latpp) + 2

# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# resolution = 'i' means use intermediate resolution coastlines.
# lon_0, lat_0 are the central longitude and latitude of the projection.
m = Basemap(
    llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax,
            resolution='i', projection='merc', lon_0=np.mean((lonmin, lonmax)), lat_0=np.mean((latmin, latmax)))
m.drawcoastlines()
m.drawcountries()

# draw parallels and meridians.
m.drawparallels(np.arange(-40, 80., 5.), color='gray')
m.drawmeridians(np.arange(-30., 80., 5.), color='gray')

# plot pierce points
print(lonpp, latpp)
x1, y1 = m(lonpp, latpp)
m.scatter(x1, y1, s=30, marker='o', color='k', alpha=.3)

plt.title('Pierce points for ' + phase + ' at ' + depth + ' km depth')
plt.show()
