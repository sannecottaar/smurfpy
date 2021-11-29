###########################
# Plot map with piercepoints
############################

#----------------------------------------#
import sys
import mpl_toolkits
import mpl_toolkits.basemap
from mpl_toolkits.basemap import Basemap
import numpy as np
import scipy
from scipy import interpolate
import matplotlib.pyplot as plt
import subprocess
import glob, sys
import pandas as pd
#----------------------------------------#

# Command line help
if len(sys.argv) != 4 or str(sys.argv[1]).lower() == 'help':
    print('\n')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print(sys.argv[0])
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('Description:           Plots discontinuity depth pierce points for 410, 660 and station locations on map.')
    print('Inputs:                path to data, filter band')
    print('Outputs:               matplotlib plot)\n')
    print('Usage:                 >> python3 plot_map_stations_pierce_points.py pathtodata rffilter')
    print('Format                 1-3: [str]')
    print('Recommended:           >> python3 plot_map_stations_pierce_points.py 410 P410s jgf1')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('\n')
    sys.exit()

def plot_pierce_points(Data, filt):
# Initial options
    rffilter = str(filt)  # RF filter
    Results = Data+'_Results'

    depths = [410, 660]
    colours = ['mediumvioletred','coral']
    shapes = ['o','s']
    plt.figure(figsize=(6, 8))

    piercelist410 = [Results+'/PP_410km_P410s_'+rffilter+'.txt']
    piercelist660 = [Results+'/PP_660km_P660s_'+rffilter+'.txt']

    # Read in pierce points
    lonpp410 = []
    latpp410 = []
    lonpp660 = []
    latpp660 = []

    for filename in piercelist410:
        print(filename)
        rd = open(filename, 'r')

        for line in rd.readlines():
            val = line.split()
            lonpp410.append(float(val[2]))
            latpp410.append(float(val[1]))

        rd.close()

    for filename in piercelist660:
        print(filename)
        rd = open(filename, 'r')

        for line in rd.readlines():
            val = line.split()
            lonpp660.append(float(val[2]))
            latpp660.append(float(val[1]))

        rd.close()

        # Set bounds of map based on piercepoints
    
    lonpp = np.concatenate((np.array(lonpp410), np.array(lonpp660)))
    latpp = np.concatenate((np.array(latpp410), np.array(latpp660)))

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
    x410, y410 = m(lonpp410, latpp410)
    m.scatter(x410, y410, s=30, marker=shapes[0], color=colours[0], alpha=.3, label='P410s')
    x660, y660 = m(lonpp660, latpp660)
    m.scatter(x660, y660, s=30, marker=shapes[1], color=colours[1], alpha=.3, label='P660s')


    data = pd.read_csv('/path/to/station/info', sep=" ")
    lonMY = (np.array(data['lon'])
    latMY = (np.array(data['lat'])

    x2, y2 = m(lonMY, latMY)
    m.scatter(x2, y2, s=30, marker='v', color='black', alpha=.9, label='NetworkName')

    plt.legend(frameon=False, loc = 2)

    plt.title('Pierce points for P410s/P660s at 410/660 km depth')
    plt.savefig(Results+'/Pierce_points_and_stations_'+filt+'.png')
 
    # plt.show()
             
plot_pierce_points(sys.argv[1], sys.argv[2])
