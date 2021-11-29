# common conversion point stacking

# import modules
from geographiclib.geodesic import Geodesic as geo
import numpy as np
import scipy
import scipy.ndimage
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import os
import os.path
import math
import msgpack
import msgpack_numpy as m
m.patch()
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
import numpy.ma as ma
import matplotlib.pyplot as plt

# definition of the half width of the fresnel zone
knotspacing = lambda r, vs: 1./2.*np.sqrt(((10./3.*vs)+r)**2.-r**2.) # in m for a 10s wave


def haversine(lat1, long1, lats2, longs2):
    """
    Calculate the distance between two points on earth in m
    """
    d = []

    for i in range(len(lats2)):
        lat2 = lats2[i]
        long2 = longs2[i]
        earth_radius = 6371.e3  # m
        dLat = math.radians(lat2 - lat1)
        dLong = math.radians(long2 - long1)

        a = (math.sin(dLat / 2) ** 2 +
             math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dLong / 2) ** 2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        d.append(earth_radius * c)
    return float(d[0])


def weight(distance, depth, vs, factor):
    '''
    Calculates the weight based on the distance and the fresnel zone width for a given depth
    '''
    delta = distance/(factor*knotspacing(depth*1.e3, vs)) # distance in m ~~~~~~ fresnel zone times factor~~~
    if delta > 2:
        weight = 0
    elif delta > 1:
        weight = .25*(2.-delta)**3.
    else:
        weight = .75*delta**3.-1.5*delta**2.+1.

    return weight

class VOL(dict):
    def __init__(self, *arg, **kw):
        super(VOL, self).__init__(*arg, **kw)
        self.__dict__ = self
    def __getattr__(self, name):
        return self[name]

class ccp_volume(object):
    """
       Handling large stacked volumes
    """
    def __init__(self, *arg, **kw):
        self.VOL = VOL(*arg, **kw)


#
# Load latest volume to dictionary
#
    def load_latest(self, name='Megavolume', filter='rff2', conversion='EU60', factor=1.):
        if os.path.isfile('../CCP_volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/Stack_master.PICKLE'):
            volumefile='../CCP_volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/Stack_master.PICKLE'
            print('loading', name, volumefile)
        else:
            line = open('../CCP_volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/filenames.dat', 'r').readlines()[-1]
    #        line = open('../AFR_CCP_volumes_SAFE/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/filenames.dat', 'r').readlines()[-1]    # A temporary solution....
            runnum = int(float(line.split()[0]))
            volumefile = line.split()[1]
            print('loading', name, runnum, volumefile)

        # Read in volumes
        self.VOL.update(msgpack.unpack(open(volumefile, 'rb'), use_list=False, object_hook=m.decode, encoding='utf-8'))
        # self.VOL.update(msgpack.unpack(open(volumefile, 'rb'), use_list=False, object_hook=m.decode))


        return self


#
# Plot data coverage map at predefined depth
#

    def plot_datacoverage(self,Data,depth,name='Megavolume',conversion='prem',filter='jgf1', factor=2.):

        coverage_file = open(Results+'/CCP_volumes/' + name + '_' + str(depth) + '_weights_' + conversion + '_' + str(filter) + '_' +str(int(factor))+'.txt', 'w')
        fig = plt.figure(figsize=(6,6))
        d = np.argmin(np.abs(self.VOL['grid_depth']-depth))
        slice = self.VOL['volumeweight'][:,:, d].copy()


        xx, yy = np.meshgrid(self.VOL['grid_lon'], self.VOL['grid_lat'])

        m = plt.axes(projection=ccrs.Mercator())

        m.add_feature(cfeature.COASTLINE)
        m.add_feature(cfeature.BORDERS, linestyle="--")

        gl = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlines = False
        gl.ylines = False
        gl.xlocator = mticker.FixedLocator([110,115,120])
        gl.ylocator = mticker.FixedLocator([0,5,10])

        x, y = xx, yy

        contours = [1., .5e1, 1.e1, .5e2, 1.e2, .5e3, 1.e3, .5e4, 1.e4]#[1.e0,1.e1,1.e2,1.e3,1.e4]
        im = plt.contourf(x, y, slice.T, contours, norm=LogNorm(),zorder=1,transform=ccrs.PlateCarree())
        # To write out the volumeweights
        for i in range(len(self.VOL['grid_lon'])):
            for j in range(len(self.VOL['grid_lat'])):
                coverage_file.write(str(self.VOL['grid_lon'][i])+'\t')
                coverage_file.write(str(self.VOL['grid_lat'][j])+'\t')
                coverage_file.write(str(slice.T[j,i])+'\n')
        coverage_file.close()

        dep_str=str(depth)+'km'

        d_lon=np.min(self.VOL['grid_lon'])+0.4
        d_lat=np.min(self.VOL['grid_lat'])+0.4
        x, y = d_lon, d_lat
        plt.text(x, y,dep_str,fontsize=18,ha='left',va='bottom',color='black', bbox=dict(facecolor='white', edgecolor='black', pad=3.0))

        fig.subplots_adjust(bottom=.2)
        cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
        cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cb.set_label('Sum of weights at ' + str(depth) + ' km', fontsize=18)
        cb.ax.tick_params(labelsize=16)



#
# Plot topography maps
#

    def plot_topography(self,mindepth,maxdepth,name='Megavolume',conversion='prem',filter='jgf1',factor=2.,mincoverage=10., amplitude = True, blobs=True):
        # Plots topography of maximum between mindepth and maxdepth, masking if sum of weights is beneath mincoverage.
        # If amplitude =True, it will plot the amplitude and not the depth
        # window buffer at end of depth array to not pick.
        wb=5

        plt.figure(figsize=(10, 8))
        depths = self.VOL['grid_depth']
        #print('depths are ', depths)
        val_list = [x for x in range(len(depths)) if depths[x] > mindepth and depths[x] < maxdepth]


        # thickness = np.empty((len(self.VOL['grid_lon']), len(self.VOL['grid_lat'])))
        dmap = np.empty((len(self.VOL['grid_lon']), len(self.VOL['grid_lat'])))
        # coverage = np.empty((len(self.VOL['grid_lon']), len(self.VOL['grid_lat'])))
        dsign = np.empty((len(self.VOL['grid_lon']), len(self.VOL['grid_lat'])))

        for i in range(len(self.VOL['grid_lon'])):
            for j in range(len(self.VOL['grid_lat'])):
                RF = self.VOL['volume'][i, j,:]/self.VOL['volumeweight'][i, j,:]
                plt.plot(RF)
                 # This appears to be doing 1.96x std so its ~2std.
                std = 1.96*np.sqrt(self.VOL['volumesigma'][i, j,:]/(self.VOL['volumeweight'][i, j,:]*self.VOL['volumeweight'][i, j,:]))
                maxmap = np.argmax(RF[val_list])

                if depths[val_list[maxmap]] > depths[val_list[wb-1]] and depths[val_list[maxmap]] < depths[val_list[-wb]]:
                    if self.VOL['volumeweight'][i, j, val_list[maxmap]] >= mincoverage:
                        if abs(RF[val_list[maxmap]])>std[val_list[maxmap]]:
                            dmap[i, j] = depths[val_list[maxmap]]


        dmapall = np.ravel(dmap)
        if amplitude == False:
            l = [l for l in range(len(dmapall)) if dmapall[l] > mindepth+1. and dmapall[l] < maxdepth - 1.]

            if np.median((dmapall[l]))*100 == np.NaN:
                print("NaN median")
            else:
                median=         round(np.median((dmapall[l]))*100)/100
                variance=       round(np.var((dmapall[l]))*100)/100

                med_str='med = ' + str(median)
                var_str='var = ' + str(variance)

                print(med_str)
                print(var_str)



        # Prepare map
        m = plt.axes(projection=ccrs.Mercator())

        
        xx, yy = np.meshgrid(self.VOL['grid_lon'], self.VOL['grid_lat'])
        x, y = xx, yy

        if amplitude is False:
            cs = plt.pcolormesh(x, y, dmap.T, vmin=mindepth, vmax=maxdepth, cmap='inferno', linewidth=0, rasterized=False,transform=ccrs.PlateCarree())
            #mask area if not significant
            dsign=ma.masked_array(dsign, dsign==0)
        #            sign = plt.contourf(x, y, dsign.T, vmin=0.5, vmax=1.0, cmap=cm.Greys, linewidth=0, alpha=0.7, rasterized=False)

        else:
            cs = plt.pcolormesh(x, y, dmap.T, vmin=0.02, vmax=0.12, cmap='inferno', linewidth=0, rasterized=False,transform=ccrs.PlateCarree())
        cs.cmap.set_under([0.8, 0.8, 0.8])
        cs.cmap.set_over([0.8, 0.8, 0.8])
        #mindepth=np.argmin(dmap.T)
        #maxdepth=np.argmax(dmap.T)

        cb = plt.colorbar(cs,ticks=[mindepth, mindepth+(maxdepth-mindepth)/6, mindepth+2*(maxdepth-mindepth)/6,mindepth+3*(maxdepth-mindepth)/6,mindepth+4*(maxdepth-mindepth)/6,mindepth+5*(maxdepth-mindepth)/6,maxdepth])

        #cb.set_label('Maximum map between ' + str(mindepth)+' and ' + str(maxdepth)+' (km)', size=30)
        if mindepth<400:
            cb.set_label('Depth of 410 (km)', size=18)
        else:
            cb.set_label('Depth of 660 (km)', size=18)
        cb.ax.tick_params(labelsize=16)
        #cb.set_label('Amplitude')

        # Label region and conversion at top.
        r_lon=np.max(self.VOL['grid_lon'])-0.4
        r_lat=np.max(self.VOL['grid_lat'])-0.4
        x, y = r_lon, r_lat
        plt.text(x, y,'region',fontsize=18,ha='right',va='top',color='black', bbox=dict(facecolor='white', edgecolor='black', pad=3.0))

        c_lon=np.min(self.VOL['grid_lon'])+0.4
        c_lat=np.max(self.VOL['grid_lat'])-0.4
        x, y = c_lon, c_lat
        plt.text(x, y,conversion,fontsize=18,ha='left',va='top',color='black', bbox=dict(facecolor='white', edgecolor='black', pad=3.0))

        # Label median and variance at bottom.
        med_lon=np.min(self.VOL['grid_lon'])+0.4
        med_lat=np.min(self.VOL['grid_lat'])+0.4
        x, y = med_lon, med_lat
        plt.text(x, y,med_str,fontsize=18,ha='left',va='bottom',color='black', bbox=dict(facecolor='white', edgecolor='black', pad=3.0))

        var_lon=np.max(self.VOL['grid_lon'])-0.4
        var_lat=np.min(self.VOL['grid_lat'])+0.4
        x, y = var_lon, var_lat
        plt.text(x, y,var_str,fontsize=18,ha='right',va='bottom',color='black', bbox=dict(facecolor='white', edgecolor='black', pad=3.0))


        # cb.set_ticks([380,400,420,440])
        cb.solids.set_rasterized(True)
        xt, yt = -13.2, 70.6

        m.add_feature(cfeature.COASTLINE)
        m.add_feature(cfeature.BORDERS, linestyle="--")

        gl = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlines = False
        gl.ylines = False
        gl.xlocator = mticker.FixedLocator([110,115,120])
        gl.ylocator = mticker.FixedLocator([0,5,10])

#
# Plot map of MTZ width
#
    def plot_mtzwidth(self,name,conversion,filter,factor, mincoverage):
        # Plots topography of maximum between mindepth and maxdepth, masking if sum of weights is beneath mincoverage.
        # window buffer at end of depth array to not pick.
        wb=5
        plt.figure(figsize=(10, 8))
        depths = self.VOL['grid_depth']
        print(depths)
        l410 = [x for x in range(len(depths)) if depths[x] > 380 and depths[x] < 430]
        l660 = [x for x in range(len(depths)) if depths[x] > 630 and depths[x] < 700]

        # thickness = np.empty((len(self.VOL['grid_lon']), len(self.VOL['grid_lat'])))
        dmap = np.empty((len(self.VOL['grid_lon']), len(self.VOL['grid_lat'])))
        # coverage = np.empty((len(self.VOL['grid_lon']), len(self.VOL['grid_lat'])))

        for i in range(len(self.VOL['grid_lon'])):
            for j in range(len(self.VOL['grid_lat'])):

                RF = self.VOL['volume'][i, j,:]/self.VOL['volumeweight'][i, j,:]
                plt.plot(RF)
                # Again originally doing 1.96x std so changes to 1xstd
                std = 1.96*np.sqrt(self.VOL['volumesigma'][i, j,:]/(self.VOL['volumeweight'][i, j,:]*self.VOL['volumeweight'][i, j,:]))
                maxmap410 = np.argmax(RF[l410])
                maxmap660 = np.argmax(RF[l660])

                if depths[l410[maxmap410]] > depths[l410[wb-1]] and depths[l410[maxmap410]] < depths[l410[-wb]] and depths[l660[maxmap660]] > depths[l660[wb-1]] and depths[l660[maxmap660]] < depths[l660[-wb]]:
                    if self.VOL['volumeweight'][i, j, l410[maxmap410]] >= mincoverage:
                        if RF[l410[maxmap410]] > std[l410[maxmap410]] and RF[l660[maxmap660]] > std[l660[maxmap660]]:
                            dmap[i, j] = depths[l660[maxmap660]]-depths[l410[maxmap410]]

        dmapall = np.ravel(dmap)

        l = [l for l in range(len(dmapall)) if dmapall[l] > 210+1. and dmapall[l] < 290 - 1.]

        median=         round(np.median((dmapall[l]))*100)/100
        variance=       round(np.var((dmapall[l]))*100)/100

        med_str='med = ' + str(median)
        var_str='var = ' + str(variance)

        print(med_str)
        print(var_str)

        # Prepare map
        m = plt.axes(projection=ccrs.Mercator())


        xx, yy = np.meshgrid(self.VOL['grid_lon'], self.VOL['grid_lat'])
        x, y = xx, yy

        #cs = plt.contourf(x, y, dmap.T, vmin=220.,levels=np.linspace(220., 280., 81.), cmap=cm.RdBu)
        cs = plt.pcolormesh(x, y, dmap.T, vmin=220., vmax=280., cmap=cm.RdBu, linewidth=0, rasterized=True,transform=ccrs.PlateCarree())


        cs.cmap.set_under([0.8, 0.8, 0.8])
        cs.cmap.set_over([0.8, 0.8, 0.8])
        cb = plt.colorbar()
        cb.set_label('Transition zone thickness (km)', size=18)
        cb.set_ticks([220,235,250,265,280])
        cb.ax.tick_params(labelsize=16)
        cb.solids.set_rasterized(True)

        dmapall = np.ravel(dmap)


        m.add_feature(cfeature.COASTLINE)
        m.add_feature(cfeature.BORDERS, linestyle="--")

        gl = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlines = False
        gl.ylines = False
        # replace with your tick locations of choice
        gl.xlocator = mticker.FixedLocator([110,115,120])
        gl.ylocator = mticker.FixedLocator([0,5,10])




    def plot_mtzwidth_write(self,name='Megavolume',conversion='prem',filter='jgf1',factor=2., mincoverage=10.):

        #This routine is also used to make a txt file with the significant depths of the 410 and the 660
        depth410660_2SE = open('../CCP_volumes/' + name + '_' + conversion+ '_'+ str(filter)+'_'+str(int(factor))+ '_MTZ_2SE_depths.txt', 'w')

        wb=5

        # plt.figure(figsize=(18, 8))
        depths = self.VOL['grid_depth']
        l410 = [x for x in range(len(depths)) if depths[x] > 380 and depths[x] < 430]
        l660 = [x for x in range(len(depths)) if depths[x] > 630 and depths[x] < 700]

        for i in range(len(self.VOL['grid_lon'])):
            for j in range(len(self.VOL['grid_lat'])):


                RF = self.VOL['volume'][i, j,:]/self.VOL['volumeweight'][i, j,:]
                # Again originally 1.96 so changed to 1.
                std = 1.96*np.sqrt(self.VOL['volumesigma'][i, j,:]/(self.VOL['volumeweight'][i, j,:]*self.VOL['volumeweight'][i, j,:]))
                max410 = np.argmax(RF[l410])
                max660 = np.argmax(RF[l660])
                maxamp410 = np.max(RF[l410])
                maxamp660 = np.max(RF[l660])

                # Store thickness, pick significance
                if RF[l410[max410]] > std[l410[max410]]:
                    sig_410=1
                else:
                    sig_410=0

                if RF[l660[max660]] > std[l660[max660]]:
                    sig_660=1
                else:
                    sig_660=0


                COV410 = False; STD410 = False; WITHINLIM410 = False;  COV660 = False; STD660 = False; WITHINLIM660 = False

                if self.VOL['volumeweight'][i,j,l410[max410]]>=mincoverage: COV410 = True
                if self.VOL['volumeweight'][i,j,l660[max660]]>=mincoverage: COV660 = True
                if RF[l410[max410]] > std[l410[max410]]                : STD410 = True
                if RF[l660[max660]] > std[l660[max660]]                : STD660 = True
                if depths[l410[max410]] > depths[l410[wb-1]] and depths[l410[max410]] < depths[l410[-wb]]: WITHINLIM410 = True
                if depths[l660[max660]] > depths[l660[wb-1]] and depths[l660[max660]] < depths[l660[-wb]]: WITHINLIM660 = True

                if COV410 and WITHINLIM410 and COV660 and WITHINLIM660:
                    if STD410 and STD660:
                        depth410660_2SE.write(str(self.VOL['grid_lat'][j])+'\t')
                        depth410660_2SE.write(str(self.VOL['grid_lon'][i])+'\t')
                        depth410660_2SE.write(str(depths[l410[max410]])+'\t')
                        depth410660_2SE.write(str(depths[l660[max660]])+'\t')
                        depth410660_2SE.write(str(sig_410)+'\t')
                        depth410660_2SE.write(str(sig_660)+'\t')
                        depth410660_2SE.write(str(maxamp410)+'\t')
                        depth410660_2SE.write(str(maxamp660)+'\n')

                    elif STD660 and not STD410:
                        depth410660_2SE.write(str(self.VOL['grid_lat'][j])+'\t')
                        depth410660_2SE.write(str(self.VOL['grid_lon'][i])+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(depths[l660[max660]])+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(sig_660)+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(maxamp660)+'\n')

                    elif STD410 and not STD660:
                        depth410660_2SE.write(str(self.VOL['grid_lat'][j])+'\t')
                        depth410660_2SE.write(str(self.VOL['grid_lon'][i])+'\t')
                        depth410660_2SE.write(str(depths[l410[max410]])+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(sig_410)+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(maxamp410)+'\t')
                        depth410660_2SE.write(str('NaN')+'\n')

                if not COV410 and COV660 and WITHINLIM660:
                    if STD660:
                        depth410660_2SE.write(str(self.VOL['grid_lat'][j])+'\t')
                        depth410660_2SE.write(str(self.VOL['grid_lon'][i])+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(depths[l660[max660]])+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(sig_660)+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(maxamp660)+'\n')

                if not COV660 and COV410 and WITHINLIM410:
                    if STD410:
                        depth410660_2SE.write(str(self.VOL['grid_lat'][j])+'\t')
                        depth410660_2SE.write(str(self.VOL['grid_lon'][i])+'\t')
                        depth410660_2SE.write(str(depths[l410[max410]])+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(sig_410)+'\t')
                        depth410660_2SE.write(str('NaN')+'\t')
                        depth410660_2SE.write(str(maxamp410)+'\t')
                        depth410660_2SE.write(str('NaN')+'\n')
                        
        #This routine is also used to make a txt file with the depth of the 410 and the 660
        depth410660_2SE.close()


#
# Plot crossections
#
    def plot_crosssection(self, direction, lonorlat, amplify=1., name='Megavolume',conversion='prem',filter='jgf1', factor=2., zoom=False, mincoverage=10):
        # set volume lats and lons
        # window buffer at end of depth array to not pick.
        wb=5

        if direction == 'NS':
            lon = lonorlat
            n = np.argmin(np.abs(self.VOL['grid_lon']-lon))
            crossec = self.VOL['volume'][n,:,:].T.copy()
            vol_sig = self.VOL['volumesigma'][n,:,:].T.copy()
            w = self.VOL['volumeweight'][n,:,:].T

            xaxis = self.VOL['grid_lat']
            xlabel = 'latitude (dg)'
            yends = [lon, lon]
            xends = [self.VOL['latmin'], self.VOL['latmax']]


        if direction == 'EW':
            lat = lonorlat
            n = np.argmin(np.abs(self.VOL['grid_lat']-lat))
            crossec = self.VOL['volume'][:, n,:].T.copy()
            vol_sig = self.VOL['volumesigma'][:, n,:].T.copy()
            w = self.VOL['volumeweight'][:, n,:].T
            xaxis = self.VOL['grid_lon']
            xlabel = 'longitude (dg)'
            xends = [self.VOL['lonmin'], self.VOL['lonmax']]
            yends = [lat, lat]

        depths = self.VOL['grid_depth']

        # normalize
        for i in range(np.shape(w)[0]):
            for j in range(np.shape(w)[1]):
                if w[i, j] > mincoverage:

                    crossec[i, j] = crossec[i, j]/w[i, j]
                    if crossec[i, j] > 0:
                        vol_sig[i, j] = crossec[i, j]-1.96*np.sqrt(vol_sig[i, j]/(w[i, j]*w[i, j]))
                        if vol_sig[i, j] < 0:
                            vol_sig[i, j] = 0.
                    if crossec[i, j] < 0:
                        vol_sig[i, j] = crossec[i, j]+1.96*np.sqrt(vol_sig[i, j]/(w[i, j]*w[i, j]))
                        if vol_sig[i, j] > 0:
                            vol_sig[i, j] = 0.

                else:
                    crossec[i, j] = 1000.



        plt.figure(figsize=(14, 8))
        plt.tight_layout()

        plt.subplot(2, 2, 2)
    
        m = plt.axes(projection=ccrs.Mercator())

        m.add_feature(cfeature.COASTLINE)
        m.add_feature(cfeature.BORDERS, linestyle="--")

        gl = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlines = False
        gl.ylines = False
        gl.xlocator = mticker.FixedLocator([110,115,120])
        gl.ylocator = mticker.FixedLocator([0,5,10])

        if direction == 'NS':
            x1, y1 = yends[0], xends[0]
            x2, y2 = yends[1], xends[1]
            m.plot([x1, x2], [y1, y2], color='r', linewidth=1, zorder=1,transform=ccrs.PlateCarree())
            x3, y3 = lon*np.ones(len(xaxis),), np.round(xaxis/10.)*10.
            m.scatter(x3, y3, 80, xaxis, zorder=2,transform=ccrs.PlateCarree())
        if direction == 'EW':
            x1, y1 = xends[0], yends[0]
            x2, y2 = xends[1], yends[1]
            m.plot([x1, x2], [y1, y2], color='r', linewidth=1, zorder=1,transform=ccrs.PlateCarree())
            x3, y3 = np.round(xaxis/10.)*10, lat*np.ones(len(xaxis),)
            m.scatter(x3, y3, 80, xaxis, zorder=2,transform=ccrs.PlateCarree())

        norm = 0.2/amplify

        # plot
        m.add_feature(cfeature.COASTLINE)
        m.add_feature(cfeature.BORDERS, linestyle="--")

        gl = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlines = False
        gl.ylines = False
        gl.xlocator = mticker.FixedLocator([110,115,120])
        gl.ylocator = mticker.FixedLocator([0,5,10])


        # plot

        plt.subplot(2, 2, 1)

        xx, yy = np.meshgrid(xaxis, depths)
        # print(xx)

        cs = plt.pcolor(xx, yy, crossec, vmin=-0.3, vmax=0.3, rasterized=True,cmap=cm.coolwarm)

        plt.colorbar()
        cs.cmap.set_over([0.8, 0.8, 0.8])
        if zoom:
            plt.ylim([300, 800])
        else:
            plt.ylim([min(depths), max(depths)])
        plt.gca().invert_yaxis()
        plt.xlim(xends)
        plt.plot([min(xaxis), max(xaxis)], [410, 410], '--k', linewidth=0.2)
        plt.plot([min(xaxis), max(xaxis)], [660, 660], '--k', linewidth=0.2)
        plt.ylabel('Depth (km)', fontsize=12)
        plt.xlabel(xlabel, fontsize=12)

        plt.subplot(2, 1, 2)


        for t in np.arange(0, len(xaxis), 1):

            lx = [x for x in range(0, len(depths)) if (np.abs(w[x, t]) > mincoverage)]# and np.abs(vol_sig[x,t])>std[x,t]/1.96)]
            RF = vol_sig[lx, t]/norm+xaxis[t]
            RFfull = crossec[lx, t]/norm+xaxis[t]
            std = 1.96*np.sqrt(vol_sig[lx, t]/(w[lx, t]*w[lx, t]))

            plt.fill_betweenx(depths[lx], RFfull, xaxis[t], where=RFfull >= xaxis[t], facecolor='k', rasterized=False)
            plt.fill_betweenx(depths[lx], RFfull, xaxis[t], where=xaxis[t] >= RFfull, facecolor='k', rasterized=False)
            plt.fill_betweenx(depths[lx], RF, xaxis[t], where=RF >= xaxis[t], facecolor=[1.0, 0., 0.], rasterized=False)
            plt.fill_betweenx(depths[lx], RF, xaxis[t], where=xaxis[t] >= RF, facecolor=[0.0, 0.0, 1.], rasterized=False)
            RF2 = crossec[lx, t]/norm
            RF3 = vol_sig[lx, t]/norm
            l410 = [x for x in range(len(depths[lx])) if depths[lx[x]] > 370 and depths[lx[x]] < 480]
            l660 = [x for x in range(len(depths[lx])) if depths[lx[x]] > 620 and depths[lx[x]] < 740]
            # This is the section where the picking occurs

            if len(l410) > 20:
                max410 = np.argmax(RF3[l410])
                ind = lx[l410[max410]]
                if depths[l410[max410]] > depths[l410[wb-1]] and depths[l410[max410]] < depths[l410[-wb]]:
                    if np.abs(RF3[l410[max410]])>std[l410[max410]]/norm: # /1.96:
                        plt.plot([xaxis[t]+0.1, 0.5*RF2[l410[max410]]+xaxis[t]], [depths[ind], depths[ind]], 'y', linewidth=4)

            if len(l660) > 20:
                max660 = np.argmax(RF3[l660])
                ind = lx[l660[max660]]
                if depths[l660[max660]] > depths[l660[wb-1]] and depths[l660[max660]] < depths[l660[-wb]]:
                    if np.abs(RF3[l660[max660]])>std[l660[max660]]/norm: # /1.96:
                        plt.plot([xaxis[t]+0.1, 0.5*RF2[l660[max660]]+xaxis[t]], [depths[ind], depths[ind]], 'y', linewidth=4)

        plt.scatter(np.round(xaxis/10.)*10., 80.*np.ones(len(xaxis),), 80, xaxis, rasterized=False)
        plt.plot([min(xaxis), max(xaxis)], [410, 410], '--k', linewidth=1)
        plt.plot([min(xaxis), max(xaxis)], [660, 660], '--k', linewidth=1)
        plt.ylabel('Depth (km)', fontsize=12)
        plt.xlabel(xlabel, fontsize=12)
        plt.xlim([min(xaxis), max(xaxis)])




        if zoom:
            plt.ylim([300, 800])
        else:
            plt.ylim([min(depths), max(depths)])
        plt.gca().invert_yaxis()
        print(min(depths))



#
# Plot crossections
#

    def plot_crosssection_any(self, lon1,lon2,lat1,lat2,numpoints=200,amplify=1.,name='Megavolume',conversion='prem',filter='jgf1',factor=2.,zoom=False,mincoverage=10.):
        # set volume lats and lons
        # window buffer at end of depth array to not pick.
        wb=5

        plt.figure(figsize=(12,8))
        plt.tight_layout()

        inv = geo.WGS84.Inverse(lat1, lon1, lat2, lon2)
        points = np.linspace(0, inv['s12'], numpoints)
        line = geo.WGS84.Line(lat1, lon1, inv['azi1'])

        lats = []
        lons = []
        for i in range(len(points)):
            lats.append(line.Position(points[i])['lat2'])
            lons.append(line.Position(points[i])['lon2'])

        lats = np.array(lats)
        lons = np.array(lons)

        crossec = []
        vol_sig = []
        w = []

        dist = []
        for i in range(len(lats)):
            dist.append(haversine(lats[0], lons[0], [lats[i]], [lons[i]])/111194.)

        # pixelize lon and lat
        row = (lons-np.min(self.VOL['grid_lon']))/(self.VOL['grid_lon'][1]-self.VOL['grid_lon'][0])
        for i in range(len(row)):
            if row[i] < 0:
                row[i] = row[i]+len(lons)
        col = (lats-np.min(self.VOL['grid_lat']))/(self.VOL['grid_lat'][1]-self.VOL['grid_lat'][0])

        for dp in range(len(self.VOL['grid_depth'])):
            crossec.append(scipy.ndimage.map_coordinates(self.VOL['volume'][:,:, dp], np.vstack((row, col))))
            vol_sig.append(scipy.ndimage.map_coordinates(self.VOL['volumesigma'][:,:, dp], np.vstack((row, col))))
            w.append(scipy.ndimage.map_coordinates(self.VOL['volumeweight'][:,:, dp], np.vstack((row, col))))


        crossec = np.array(crossec)
        vol_sig = np.array(vol_sig)
        w = np.array(w)

        xaxis = self.VOL['grid_lat']
        xends = [lon1, lon2]
        yends = [lat1, lat2]

        depths = self.VOL['grid_depth']

        # normalize
        for i in range(np.shape(w)[0]):
            for j in range(np.shape(w)[1]):
                if w[i, j] > mincoverage:

                    crossec[i, j] = crossec[i, j]/w[i, j]
                    if crossec[i, j] > 0:
                        vol_sig[i, j] = crossec[i, j]-1.96*np.sqrt(vol_sig[i, j]/(w[i, j]*w[i, j]))
                        if vol_sig[i, j] < 0:
                            vol_sig[i, j] = 0.
                    if crossec[i, j] < 0:
                        vol_sig[i, j] = crossec[i, j]+1.96*np.sqrt(vol_sig[i, j]/(w[i, j]*w[i, j]))
                        if vol_sig[i, j] > 0:
                            vol_sig[i, j] = 0.
                else:
                    crossec[i, j] = 1000.

        # Map
        plt.subplot(2, 2, 2)
        
        m = plt.axes(projection=ccrs.Mercator())

        x1, y1 = xends[0], yends[0]
        x2, y2 = xends[1], yends[1]
        m.plot([x1, x2], [y1, y2], color='r', linewidth=1, zorder=1,transform=ccrs.PlateCarree())

        m.add_feature(cfeature.COASTLINE)
        m.add_feature(cfeature.BORDERS, linestyle="--")

        gl = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlines = False
        gl.ylines = False
        gl.xlocator = mticker.FixedLocator([110,115,120])
        gl.ylocator = mticker.FixedLocator([0,5,10])


        # Stack 1: colourmap
        plt.subplot(2, 2, 1)
        xx, yy = np.meshgrid(dist, depths)

        cs = plt.pcolormesh(xx, yy, crossec, vmin=-0.3, vmax=0.3, rasterized=True,cmap=cm.coolwarm)
        plt.colorbar()
        cs.cmap.set_over([0.8, 0.8, 0.8])
        if zoom:
            plt.ylim([300, 800])
        else:
            plt.ylim([min(depths), 400])
        plt.gca().invert_yaxis()
        plt.xlim([min(dist),max(dist)])
        plt.plot([min(xaxis), max(xaxis)], [410, 410], '--k', linewidth=0.2)
        plt.plot([min(xaxis), max(xaxis)], [660, 660], '--k', linewidth=0.2)

        # corrected by 3D model
        # normalize

        norm = 0.2/amplify#np.max(np.max(np.abs(crossec_3D)))/amplify

        plt.ylabel('Depth (km)', fontsize=18)
        plt.xlabel('Distance', fontsize=18)


        ax=plt.subplot(2,1,2)
        pos1 = ax.get_position() # get the original position
        pos2 = [pos1.x0 , pos1.y0 ,  pos1.width/14*8, pos1.height]
        pos2 = [pos1.x0 , pos1.y0 ,  pos1.width, pos1.height]
        ax.set_position(pos2) # set a new position

        # plot

        for t in np.arange(0, len(dist), 1):

            lx = [x for x in range(len(depths)) if (np.abs(w[x, t]) > mincoverage)]# and np.abs(vol_sig[x,t])>std[x,t]/1.96)]
            RF = vol_sig[lx, t]/norm+dist[t]
            RFfull = crossec[lx, t]/norm+dist[t]
            std = 1.96*np.sqrt(vol_sig[lx, t]/(w[lx, t]*w[lx, t]))

            plt.fill_betweenx(depths[lx], RFfull, dist[t], where=RFfull >= dist[t], facecolor='k',           rasterized=False)
            plt.fill_betweenx(depths[lx], RFfull, dist[t], where=dist[t] >= RFfull, facecolor='k',           rasterized=False)
            plt.fill_betweenx(depths[lx], RF,     dist[t], where=RF >= dist[t],     facecolor=[1.0, 0., 0.], rasterized=False)
            plt.fill_betweenx(depths[lx], RF,     dist[t], where=dist[t] >= RF,     facecolor=[0.0, 0.0, 1.],rasterized=False)
            RF2 = crossec[lx, t]/norm
            RF3 = vol_sig[lx, t]/norm
            l410 = [x for x in range(len(depths[lx])) if depths[lx[x]] > 370 and depths[lx[x]] < 480]
            l660 = [x for x in range(len(depths[lx])) if depths[lx[x]] > 620 and depths[lx[x]] < 740]
            # This is the section where the picking occurs
            # Need to add the std and array end catches. Min coverage is dealt with above.

            if len(l410) > 20:
                max410 = np.argmax(RF3[l410])
                ind = lx[l410[max410]]
                if depths[l410[max410]] > depths[l410[wb-1]] and depths[l410[max410]] < depths[l410[-wb]]:
                    if np.abs(RF3[l410[max410]])>std[l410[max410]]/norm: # /1.96:
                        plt.plot([dist[t]+0.1, 0.5*RF2[l410[max410]]+dist[t]], [depths[ind], depths[ind]], 'y', linewidth=4)

            if len(l660) > 20:
                max660 = np.argmax(RF3[l660])
                ind = lx[l660[max660]]
                if depths[l660[max660]] > depths[l660[wb-1]] and depths[l660[max660]] < depths[l660[-wb]]:
                    if np.abs(RF3[l660[max660]])>std[l660[max660]]/norm: # /1.96:
                        plt.plot([dist[t]+0.1, 0.5*RF2[l660[max660]]+dist[t]], [depths[ind], depths[ind]], 'y', linewidth=4)

        plt.ylabel('Depth (km)', fontsize=12)
        plt.xlabel('Angular distance (dg)', fontsize=12)
        plt.xlim([min(dist), max(dist)])
        plt.plot([min(dist), max(dist)], [410, 410], '--k', linewidth=1)
        plt.plot([min(dist), max(dist)], [660, 660], '--k', linewidth=1)
        if zoom:
            plt.ylim([300, 800])
        else:
            plt.ylim([min(depths), max(depths)])
        plt.gca().invert_yaxis()



#
# Plot moveout aligned on 410 or 660
#

    def plot_moveout(self,d660=True,name='Megavolume',conversion='prem',filter='jgf1',factor=2.):
    # Picks all profiles in grid and organizes them by 660 depth
        plt.figure(figsize=(8, 8))
        depths = self.VOL['grid_depth']

        if d660:
            ldisc = [x for x in range(len(depths)) if depths[x] > 620 and depths[x] < 700]
            disc = np.arange(620., 700., 2.)
            ldisc2 = [x for x in range(len(depths)) if depths[x] > 370 and depths[x] < 450]
        else:
            ldisc = [x for x in range(len(depths)) if depths[x] > 370 and depths[x] < 450]
            disc = np.arange(370, 450, 2.)
            ldisc2 = [x for x in range(len(depths)) if depths[x] > 620 and depths[x] < 700]




        weights = np.zeros((len(disc),))
        moveout = np.zeros((len(disc), len(depths)))
        d660l = []
        d410l = []
        for i in range(len(self.VOL['grid_lon'])):
            for j in range(len(self.VOL['grid_lat'])):
                RF = self.VOL['volume'][i, j,:].copy()
                for k in range(len(depths)):
                    if self.VOL['volumeweight'][i, j, k] > 0:
                        RF[k] = RF[k]/self.VOL['volumeweight'][i, j, k]

                std = 1.96*np.sqrt(self.VOL['volumesigma'][i, j,:]/(self.VOL['volumeweight'][i, j,:]*self.VOL['volumeweight'][i, j,:]))

                maxdisc = np.argmax(RF[ldisc])


                if RF[ldisc[maxdisc]] > std[ldisc[maxdisc]] :
                    d660l.append(depths[ldisc[maxdisc]])
                    n = np.nanargmin(np.abs(depths[ldisc[maxdisc]]-disc))
                    weights[n] = weights[n]+1.
                    moveout[n,:] = moveout[n,:]+RF
#                    if d660:
                    maxdisc = np.nanargmax(RF[ldisc2])
                    d410l.append(depths[ldisc2[maxdisc]])


        d660l = np.array(d660l)
        d410l = np.array(d410l)

        d660c = np.arange(620., 700., 2.)
        d410c = []
        err = []
        for i in range(len(d660c)):
            filt = [d410l[x] for x in range(len(d410l)) if (d660l[x] == d660c[i] and d410l[x] > 370 and d410l[x] < 450)]
            d410c.extend([np.mean(filt)])
            err.extend([np.std(filt)])
        d410c = np.array(d410c)


        midp = []
        discp = []
        discp2 = []
        for i in range(len(disc)):
            if weights[i] > 0:
                moveout[i,:] = moveout[i,:]/float(weights[i])

                if d660:
                    max660 = np.argmax(moveout[i, ldisc])
                    discp.append(depths[ldisc[max660]])
                    midp.append(disc[i])
                    l410 = [x for x in range(len(depths)) if depths[x] > 370 and depths[x] < 450]
                    max410 = np.argmax(moveout[i, l410])
                    discp2.append(depths[l410[max410]])
                else:
                    max410 = np.argmax(moveout[i, ldisc])
                    discp.append(depths[ldisc[max410]])
                    midp.append(disc[i])
                    l660 = [x for x in range(len(depths)) if depths[x] > 620 and depths[x] < 700]
                    max660 = np.argmax(moveout[i, l660])
                    discp2.append(depths[l660[max660]])


        ax = plt.subplot2grid((3, 6), (0, 0), colspan=5)

        if d660:
            n, bins, patches = plt.hist(d660l, np.arange(620., 700., 2.), histtype='bar')
            cmbl = plt.cm.get_cmap('bone_r')
            for d, p in zip(n, patches):
                plt.setp(p, 'facecolor', cmbl(d/max(n)*0.6+0.2))
            plt.text(646, 180, 'a.', fontsize=12)
            plt.gca().set_xticks([650, 660, 670, 680, 690])
            plt.xlim([645, 695])
        else:
            n, bins, patches = plt.hist(d660l, np.arange(370., 450., 2.), histtype='bar')
            cmbl = plt.cm.get_cmap('bone_r')
            for d, p in zip(n, patches):
                plt.setp(p, 'facecolor', cmbl(d/max(n)*0.6+0.2))
            plt.text(386, 180, 'a.', fontsize=12)
            plt.gca().set_xticks([390, 400, 410, 420, 430])
            plt.xlim([385, 435])

        plt.ylabel('# of grid points')
        plt.gca().set_yticks([0, 50, 100, 150, 200])
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax = plt.subplot2grid((3, 6), (1, 0), rowspan=2, colspan=5)

        cs = plt.pcolor(disc-1., depths, moveout.T, vmin=-0.1, vmax=0.1, rasterized=True, cmap=cm.RdYlBu_r)

        plt.plot(midp, discp, '--k', linewidth=2)

        plt.plot(midp, discp2, '--k', linewidth=2)

        plt.ylabel('Depth (km)')

        if d660:
            plt.xlim([645, 695])
            plt.xlabel('660 depth (km)')
            plt.text(646, 380, 'b.', fontsize=12)
        else:
            plt.xlim([385, 435])
            plt.xlabel('410 depth (km)')
            plt.text(386, 380, 'b.', fontsize=12)
        plt.ylim([350, 750])
        plt.gca().invert_yaxis()
        box = plt.gca().get_position()
        axcb = plt.axes([box.x0*1.05+box.width*1.05, box.y0+0.1, 0.01, box.height-0.2])
        cb = plt.colorbar(cs, cax=axcb, orientation='vertical')
        cb.set_ticks([-0.1, 0., .1])
        cb.set_label('relative amplitude')

        plt.suptitle('region' + ' - ' + str(conversion))



def plot_rectangle(bmap, lonmin,lonmax,latmin,latmax, col):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    bmap.plot(xs, ys,latlon = True, color=col)

def plot_line (bmap,lonmin,lonmax,latmin,latmax,col):
    xs = [lonmax,lonmax,lonmin]
    ys = [latmin,latmin,latmax]
    bmap.plot(xs, ys,latlon = True, color="Black")
