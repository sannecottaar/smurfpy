# common conversion point stacking

# import modules
import sys
from geographiclib.geodesic import Geodesic as geo
import numpy as np
# from matplotlib.mlab import griddata
import scipy
import scipy.ndimage
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os
import os.path
import math
import msgpack
import msgpack_numpy as m
m.patch()
import shutil
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from obspy import read
from scipy import stats
import numpy.ma as ma

import mpl_toolkits

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
        line = open('../CCP_volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/filenames.dat', 'r').readlines()[-1]
        runnum = int(float(line.split()[0]))
        volumefile = line.split()[1]
        print('loading', name, runnum, volumefile)



        # get last stackfile name

        # Read in volumes

        self.VOL.update(msgpack.unpack(open(volumefile, 'rb'), use_list=False, object_hook=m.decode, encoding='utf-8'))
        # del self.VOL['trackRFs']


        return self


#
# Plot crossections
#

    def plot_crosssection(self, direction, lonorlat, amplify=1., name='Megavolume', filter='rff2', conversion='EU60', factor=2., zoom=False, mincoverage=10):
        # set volume lats and lons

        if direction == 'NS':
            lon = lonorlat
            n = np.argmin(np.abs(self.VOL.grid_lon-lon))
            crossec = self.VOL.volume[n,:,:].T.copy()
            vol_sig = self.VOL.volumesigma[n,:,:].T.copy()
            w = self.VOL.volumeweight[n,:,:].T

            xaxis = self.VOL.grid_lat
            xlabel = 'latitude (dg)'
            yends = [lon, lon]
            xends = [self.VOL.latmin, self.VOL.latmax]


            lons = lon*np.ones_like(xaxis)
            lats = xaxis

        if direction == 'EW':
            lat = lonorlat
            n = np.argmin(np.abs(self.VOL.grid_lat-lat))
            crossec = self.VOL.volume[:, n,:].T.copy()
            vol_sig = self.VOL.volumesigma[:, n,:].T.copy()
            w = self.VOL.volumeweight[:, n,:].T
            xaxis = self.VOL.grid_lon
            xlabel = 'longitude (dg)'
            xends = [self.VOL.lonmin, self.VOL.lonmax]
            yends = [lat, lat]
            lons = xaxis
            lats = lat*np.ones_like(xaxis)

        depths = self.VOL.grid_depth

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

        plt.subplot(2, 2, 2)
        m = Basemap(projection='merc', llcrnrlat=self.VOL.latmin, urcrnrlat=self.VOL.latmax, llcrnrlon=self.VOL.lonmin, urcrnrlon=self.VOL.lonmax, lat_ts=20, resolution='i')
        m.drawparallels(np.arange(0, 70, 10.), labels=[1, 0, 0, 1], labelstyle='+/-', fontsize=10)
        m.drawmeridians(np.arange(-10, 60, 10.), labels=[1, 0, 0, 1], labelstyle='+/-', fontsize=10)

        m.drawcoastlines()
        m.drawcountries()

        m.drawmapboundary(fill_color=[1.0, 1.0, 1.0])
        if direction == 'NS':
            x1, y1 = m(yends[0], xends[0])
            x2, y2 = m(yends[1], xends[1])
            m.plot([x1, x2], [y1, y2], color='r', linewidth=1, zorder=1)
            
            x3, y3 = m(lon*np.ones(len(xaxis),), np.round(xaxis/10.)*10.)
            m.scatter(x3, y3, 80, xaxis, zorder=2)
        if direction == 'EW':
            x1, y1 = m(xends[0], yends[0])
            x2, y2 = m(xends[1], yends[1])
            m.plot([x1, x2], [y1, y2], color='r', linewidth=1, zorder=1)
            x3, y3 = m(np.round(xaxis/10.)*10, lat*np.ones(len(xaxis),)  )
            m.scatter(x3, y3, 80, xaxis, zorder=2)



        norm = 0.2/amplify

        # plot
        
        plt.subplot(2, 2, 1)

        xx, yy = np.meshgrid(xaxis, depths)
        print(xx)

        cs = plt.pcolor(xx, yy, crossec, vmin=-0.3, vmax=0.3, rasterized=True,cmap=cm.coolwarm)

        #cs = plt.pcolor(xx, yy, crossec, vmin=-0.15, vmax=0.15, rasterized=True,cmap=cm.coolwarm)
        plt.colorbar()
        cs.cmap.set_over([0.8, 0.8, 0.8])
        if zoom:
            plt.ylim([300, 800])
        else:
            plt.ylim([min(depths), max(depths)])
        plt.gca().invert_yaxis()
        plt.xlim(xends)
 
        plt.subplot(2, 1, 2)

        
        for t in np.arange(0, len(xaxis), 1):

            lx = [x for x in range(0, len(depths)) if (np.abs(w[x, t]) > mincoverage)]# and np.abs(vol_sig[x,t])>std[x,t]/1.96)]
            RF = vol_sig[lx, t]/norm+xaxis[t]
            RFfull = crossec[lx, t]/norm+xaxis[t]
            
            plt.fill_betweenx(depths[lx], RFfull, xaxis[t], where=RFfull >= xaxis[t], facecolor='k', rasterized=True)
            plt.fill_betweenx(depths[lx], RFfull, xaxis[t], where=xaxis[t] >= RFfull, facecolor='k', rasterized=True)
            plt.fill_betweenx(depths[lx], RF, xaxis[t], where=RF >= xaxis[t], facecolor=[1.0, 0., 0.], rasterized=True)
            plt.fill_betweenx(depths[lx], RF, xaxis[t], where=xaxis[t] >= RF, facecolor=[0.0, 0.0, 1.], rasterized=True)
        plt.scatter(np.round(xaxis/10.)*10., 80.*np.ones(len(xaxis),), 80, xaxis, rasterized=True)
        plt.plot([-180, 140], [410, 410], '--k', linewidth=2)
        #plt.plot([-180, 140], [520, 520], '--k', linewidth=2)
        plt.plot([-180, 140], [660, 660], '--k', linewidth=2)
        plt.ylabel('Depth (km)')
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

    def plot_crosssection_any(self, lon1,lon2,lat1,lat2,numpoints=200,amplify=1.,name='Megavolume',filter='rff2', conversion='EU60', factor=2.,zoom=False,mincoverage=10.):
        # set volume lats and lons

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
        row = (lons-np.min(self.VOL.grid_lon))/(self.VOL.grid_lon[1]-self.VOL.grid_lon[0])
        for i in range(len(row)):
            if row[i] < 0:
                row[i] = row[i]+len(self.lon)
        col = (lats-np.min(self.VOL.grid_lat))/(self.VOL.grid_lat[1]-self.VOL.grid_lat[0])

        for dp in range(len(self.VOL.grid_depth)):
            crossec.append(scipy.ndimage.map_coordinates(self.VOL.volume[:,:, dp], np.vstack((row, col))))
            vol_sig.append(scipy.ndimage.map_coordinates(self.VOL.volumesigma[:,:, dp], np.vstack((row, col))))
            w.append(scipy.ndimage.map_coordinates(self.VOL.volumeweight[:,:, dp], np.vstack((row, col))))


        crossec = np.array(crossec)
        vol_sig = np.array(vol_sig)
        w = np.array(w)

        xaxis = self.VOL.grid_lat
        xlabel = 'latitude (dg)'
        xends = [lon1, lon2]
        yends = [lat1, lat2]

        depths = self.VOL.grid_depth

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
                    crossec[i, j] = 100.
                   

        plt.subplot(2, 2, 2)
        m = Basemap(projection='merc', llcrnrlat=self.VOL.latmin, urcrnrlat=self.VOL.latmax, llcrnrlon=self.VOL.lonmin, urcrnrlon=self.VOL.lonmax, lat_ts=20, resolution='i')
        m.drawparallels(np.arange(0, 90, 10.), labels=[1, 0, 0, 1], labelstyle='+/-', fontsize=10)
        m.drawmeridians(np.arange(-180, 60, 10.), labels=[1, 0, 0, 1], labelstyle='+/-', fontsize=10)

        m.drawcoastlines()
        m.drawcountries()

        m.drawmapboundary(fill_color=[1.0, 1.0, 1.0])
        x1, y1 = m(xends[0], yends[0])
        x2, y2 = m(xends[1], yends[1])
        m.plot([x1, x2], [y1, y2], color='r', linewidth=1, zorder=1)


        plt.subplot(2, 2, 1)
        xx, yy = np.meshgrid(dist, depths)

        cs = plt.pcolor(xx, yy, crossec, vmin=-0.15, vmax=0.15, rasterized=True,cmap=cm.coolwarm)
        plt.colorbar()
        cs.cmap.set_over([0.8, 0.8, 0.8])
        if zoom:
            plt.ylim([200, 900])
        else:
            plt.ylim([min(depths), max(depths)])
        plt.gca().invert_yaxis()
  


        # corrected by 3D model
        # normalize

        norm = 0.2/amplify#np.max(np.max(np.abs(crossec_3D)))/amplify


        ax=plt.subplot(2,1,2)
        pos1 = ax.get_position() # get the original position 
        pos2 = [pos1.x0 , pos1.y0 ,  pos1.width/14*8, pos1.height]
        pos2 = [pos1.x0 , pos1.y0 ,  pos1.width, pos1.height]
        print(pos1,pos2)
        ax.set_position(pos2) # set a new position
        print(ax.get_position())#

            
        print(ax.get_position())


    # plot

        for t in np.arange(0, len(dist), 1):

            lx = [x for x in range(len(depths)) if (np.abs(w[x, t]) > mincoverage)]# and np.abs(vol_sig[x,t])>std[x,t]/1.96)]
            RF = vol_sig[lx, t]/norm+dist[t]
            RFfull = crossec[lx, t]/norm+dist[t]

            plt.fill_betweenx(depths[lx], RFfull, dist[t], where=RFfull >= dist[t], facecolor='k', rasterized=True)
            plt.fill_betweenx(depths[lx], RFfull, dist[t], where=dist[t] >= RFfull, facecolor='k', rasterized=True)
            plt.fill_betweenx(depths[lx], RF, dist[t], where=RF >= dist[t], facecolor=[1.0, 0., 0.], rasterized=True)
            plt.fill_betweenx(depths[lx], RF, dist[t], where=dist[t] >= RF, facecolor=[0.0, 0.0, 1.], rasterized=True)
            RF2 = crossec[lx, t]/norm
            l410 = [x for x in range(len(depths[lx])) if depths[lx[x]] > 366 and depths[lx[x]] < 454]
            l660 = [x for x in range(len(depths[lx])) if depths[lx[x]] > 616 and depths[lx[x]] < 704]
            if len(l410) > 20:
                max410 = np.argmax(RF2[l410])
                ind = lx[l410[max410]]
                plt.plot([dist[t]+0.1, 0.5*RF2[l410[max410]]+dist[t]], [depths[ind], depths[ind]], 'y', linewidth=2)
            if len(l660) > 20:
                max660 = np.argmax(RF2[l660])
                ind = lx[l660[max660]]
                plt.plot([dist[t]+0.1, 0.5*RF2[l660[max660]]+dist[t]], [depths[ind], depths[ind]], 'y', linewidth=2)

        plt.ylabel('Depth (km)', fontsize=28)
        plt.xlabel('Angular distance (dg)', fontsize=28)
        plt.xlim([min(dist), max(dist)])
        plt.plot([-5, 40], [410, 410], '--k', linewidth=2)
        plt.plot([-5, 40], [660, 660], '--k', linewidth=2)
        if zoom:
            plt.ylim([200, 800])
        else:
            plt.ylim([min(depths), max(depths)])
        plt.gca().invert_yaxis()
        print(ax.get_position())    



#
# Plot data coverage map at predefined depth
#

    def plot_datacoverage(self,depth,name='Megavolume',filter='rff2', conversion='EU60', factor=2.):
 
        fig = plt.figure(figsize=(6,6))
        d = np.argmin(np.abs(self.VOL.grid_depth-depth))
        slice = self.VOL.volumeweight[:,:, d].copy()


        xx, yy = np.meshgrid(self.VOL.grid_lon, self.VOL.grid_lat)

        m = Basemap(projection='merc', llcrnrlat=np.min(self.VOL.grid_lat), urcrnrlat=np.max(self.VOL.grid_lat), llcrnrlon=np.min(self.VOL.grid_lon), urcrnrlon=np.max(self.VOL.grid_lon), lat_ts=20, resolution='i')
        m.drawparallels(np.arange(0, 80, 10.), labels=[1, 0, 0, 1], linewidth=0.5, dashes=[4, 2], labelstyle='+/-', fontsize=22)
        m.drawmeridians(np.arange(-160, -130, 20.), labels=[1, 0, 0, 1], linewidth=0.5, dashes=[4, 2], labelstyle='+/-', fontsize=22)
        m.drawcountries()
        coasts = m.drawcoastlines(zorder=2, color='k', linewidth=1)

        m.drawmapboundary(fill_color=[1.0, 1.0, 1.0])
        x, y = m(xx, yy)

        contours = [1., 1.e1, 1.e2, 1.e3, 1.e4]#[1.e0,1.e1,1.e2,1.e3,1.e4]
        im =plt.contourf(x, y, slice.T, contours, norm=LogNorm(),zorder=1)

        lonmin=-157
        lonmax=-143
        latmin=58
        latmax=67
        #plot_line (m,lonmin,lonmax,latmin,latmax,col="Black")

        lonmin=-159
        lonmax=-131
        latmin=59.5
        latmax=66.5
        plot_line (m,lonmin,lonmax,latmin,latmax,col="Black")#BB

        lonmin=-160
        lonmax=-142
        latmin=58
        latmax=63
        plot_line (m,lonmin,lonmax,latmin,latmax,col="Black")#CC

        lonmin=-157
        lonmax=-140
        latmin=57
        latmax=69
        plot_line (m,lonmin,lonmax,latmin,latmax,col="Black")#AA

        lonmin=-164
        lonmax=-154
        latmin=56.5
        latmax=60.5
        plot_line (m,lonmin,lonmax,latmin,latmax,col="Black")#DD

        

        lonmin = -155.5
        lonmax = -152
        latmin =  66
        latmax =  68
        #plot_rectangle(m, lonmin,lonmax,latmin,latmax, col="Black")

        lonmin = -154
        lonmax = -148.5
        latmin =  58
        latmax =  60.5
        #plot_rectangle(m, lonmin,lonmax,latmin,latmax, col="Black")

        fig.subplots_adjust(bottom=.2)
        cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
        cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cb.set_label('Sum of weights at ' + str(depth) + ' km', fontsize=36)
    


#
# Plot topography maps
#

    def plot_topography(self,mindepth,maxdepth,name='Megavolume',filter='rff2',conversion='prem',factor=2.,mincoverage=10., amplitude = True, blobs=True):
        # Plots topography of maximum between mindepth and maxdepth, masking if sum of weights is beneath mincoverage.
        # If amplitude =True, it will plot the amplitude and not the depth
        
        plt.figure(figsize=(10, 8))
        depths = self.VOL.grid_depth
        #print('depths are ', depths)
        val_list = [x for x in range(len(depths)) if depths[x] > mindepth and depths[x] < maxdepth]
        

        thickness = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        dmap = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        coverage = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        dsign = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        
        amparray=[]
        for i in range(len(self.VOL.grid_lon)):
            for j in range(len(self.VOL.grid_lat)):

                RF = self.VOL.volume[i, j,:]/self.VOL.volumeweight[i, j,:]
                plt.plot(RF)
                
                std = 1.96*np.sqrt(self.VOL.volumesigma[i, j,:]/(self.VOL.volumeweight[i, j,:]*self.VOL.volumeweight[i, j,:]))
                maxmap = np.argmax(RF[val_list])
                
                if amplitude == False:
                    dmap[i, j] = depths[val_list[maxmap]]
                else:
                    dmap[i,j] = RF[val_list[maxmap]]

                #array for significance, 0=significant, 1=not significant    
                if abs(RF[val_list[maxmap]])>std[val_list[maxmap]]:
                    dsign[i,j] = 0.0
                else:
                    dsign[i,j] = 1.

                    
                if self.VOL.volumeweight[i, j, val_list[maxmap]] < mincoverage:
                    dmap[i, j] = 1000.
                    dsign[i,j]=0.

        #plt.show()



        # Prepare map
        m = Basemap(projection='merc', llcrnrlat=np.min(self.VOL.grid_lat), urcrnrlat=np.max(self.VOL.grid_lat), llcrnrlon=np.min(self.VOL.grid_lon), urcrnrlon=np.max(self.VOL.grid_lon), lat_ts=20, resolution='i')

        m.drawparallels(np.arange(0, 90, 5.), labels=[1, 0, 0, 1], linewidth=0.5, dashes=[4, 2], labelstyle='+/-', fontsize=20)
        m.drawmeridians(np.arange(-180, -110, 10.), labels=[1, 0, 0, 1], linewidth=0.5, dashes=[4, 2], labelstyle='+/-', fontsize=20)
        m.drawcountries()
        coasts = m.drawcoastlines(zorder=1, color='k', linewidth=1)




        xx, yy = np.meshgrid(self.VOL.grid_lon, self.VOL.grid_lat)
        x, y = m(xx, yy)

        if amplitude is False:
            cs = plt.pcolor(x, y, dmap.T, vmin=mindepth, vmax=maxdepth, cmap=cm.BrBG, linewidth=0, rasterized=True)
            #mask area if not significant
            dsign=ma.masked_array(dsign, dsign==0)
            sign = plt.contourf(x, y, dsign.T, vmin=0.5, vmax=1.0, cmap=cm.Greys, linewidth=0, alpha=0.7, rasterized=False)
        
        else:
            cs = plt.pcolor(x, y, dmap.T, vmin=0.02, vmax=0.12, cmap=cm.pink_r, linewidth=0, rasterized=True)
        cs.cmap.set_under([0.8, 0.8, 0.8])
        cs.cmap.set_over([0.8, 0.8, 0.8])
        #mindepth=np.argmin(dmap.T)
        #maxdepth=np.argmax(dmap.T)
        
        cb = plt.colorbar(cs,ticks=[mindepth, mindepth+(maxdepth-mindepth)/6, mindepth+2*(maxdepth-mindepth)/6,mindepth+3*(maxdepth-mindepth)/6,mindepth+4*(maxdepth-mindepth)/6,mindepth+5*(maxdepth-mindepth)/6,maxdepth])

        #cb.set_label('Maximum map between ' + str(mindepth)+' and ' + str(maxdepth)+' (km)', size=30)
        if mindepth==380:
            cb.set_label('Depth of 410 (km)', size=30)
        else:
            cb.set_label('Depth of 660 (km)', size=30)
        cb.ax.tick_params(labelsize=30)
        #cb.set_label('Amplitude')
        
        # cb.set_ticks([380,400,420,440])
        cb.solids.set_rasterized(True)
        xt, yt = m(-13.2, 70.6)

        m.drawcoastlines(zorder=1, color='k', linewidth=1)
        dmapall = np.ravel(dmap)
        if amplitude == False:
            l = [l for l in range(len(dmapall)) if dmapall[l] > mindepth+1. and dmapall[l] < maxdepth - 1.]
            print ('median', np.median((dmapall[l])))
            print ('variance', np.var((dmapall[l])))


        if blobs == True:
            latblob = []
            lonblob = []
            with open ('/raid1/annemijn/scripts/CCP/areared.txt') as blobs:
                for line in blobs:
                    row = line.split()
                    latblob.append(float(row[0]))
                    lonblob.append(float(row[1]))
            print (latblob)
            x,y = m(lonblob, latblob)
            m.scatter(x,y,marker='o',color='firebrick')

            latblob = []
            lonblob = []
            with open ('/raid1/annemijn/scripts/CCP/areayellow.txt') as blobs:
                for line in blobs:
                    row = line.split()
                    latblob.append(float(row[0]))
                    lonblob.append(float(row[1]))
            print (latblob)
            x,y = m(lonblob, latblob)
            m.scatter(x,y,marker='o',color='gold')

            latblob = []
            lonblob = []
            with open ('/raid1/annemijn/scripts/CCP/areagreen.txt') as blobs:
                for line in blobs:
                    row = line.split()
                    latblob.append(float(row[0]))
                    lonblob.append(float(row[1]))
            print (latblob)
            x,y = m(lonblob, latblob)
            m.scatter(x,y,marker='o',color='forestgreen')

            latblob = []
            lonblob = []
            with open ('/raid1/annemijn/scripts/CCP/areablue.txt') as blobs:
                for line in blobs:
                    row = line.split()
                    latblob.append(float(row[0]))
                    lonblob.append(float(row[1]))
            print (latblob)
            x,y = m(lonblob, latblob)
            m.scatter(x,y,marker='o',color='dodgerblue')


        lonmin = -142.1
        lonmax = -128
        latmin =  69
        latmax =  74
        #plot_rectangle(m, lonmin,lonmax,latmin,latmax)

        lonmin = -136
        lonmax = -128
        latmin =  54
        latmax =  63
        #plot_rectangle(m, lonmin,lonmax,latmin,latmax)

        lonmin = -157
        lonmax = -150
        latmin =  67.5
        latmax =  69
        #plot_rectangle(m, lonmin,lonmax,latmin,latmax)

        lonmin = -154
        lonmax = -149
        latmin =  66
        latmax =  68
        #plot_rectangle(m, lonmin,lonmax,latmin,latmax, col="Black")

        lonmin = -155
        lonmax = -150
        latmin =  56.5
        latmax =  59.5
        #plot_rectangle(m, lonmin,lonmax,latmin,latmax, col="Black")

        lonmin=-159
        lonmax=-131
        latmin=59.5
        latmax=66.5
        #plot_line (m,lonmin,lonmax,latmin,latmax,col="Black")#BB

        lonmin=-160
        lonmax=-142
        latmin=58
        latmax=63
        #plot_line (m,lonmin,lonmax,latmin,latmax,col="Black")#CC

        lonmin=-157
        lonmax=-140
        latmin=57
        latmax=69
        #plot_line (m,lonmin,lonmax,latmin,latmax,col="Black")#AA

#
# Plot map of MTZ width
#





    def plot_mtzwidth(self,name='Megavolume',filter='rff2',conversion='EU60',factor=2., mincoverage=10.):
        # Plots topography of maximum between mindepth and maxdepth, masking if sum of weights is beneath mincoverage.

        #depth410660 = open('/raid3/annemijn/scripts/CCP/mincov20MTZ_'+conversion+'_'+filter+'_'+str(int(factor))+'.txt', 'w')
    
        plt.figure(figsize=(10, 8))
        depths = self.VOL.grid_depth
        print(depths)
        l410 = [x for x in range(len(depths)) if depths[x] > 380 and depths[x] < 430]
        l660 = [x for x in range(len(depths)) if depths[x] > 630 and depths[x] < 700]

        thickness = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        dmap = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        coverage = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        
        amparray=[]
        for i in range(len(self.VOL.grid_lon)):
            for j in range(len(self.VOL.grid_lat)):

                RF = self.VOL.volume[i, j,:]/self.VOL.volumeweight[i, j,:]
                plt.plot(RF)
                
                std = 1.96*np.sqrt(self.VOL.volumesigma[i, j,:]/(self.VOL.volumeweight[i, j,:]*self.VOL.volumeweight[i, j,:]))
                maxmap410 = np.argmax(RF[l410])
                maxmap660 = np.argmax(RF[l660])
                
                if RF[l410[maxmap410]] > std[l410[maxmap410]] and RF[l660[maxmap660]] > std[l660[maxmap660]]:
                    dmap[i, j] = depths[l660[maxmap660]]-depths[l410[maxmap410]]

                #depth410660.write(str(self.VOL.grid_lat[j])+'\t')
                #depth410660.write(str(self.VOL.grid_lon[i])+'\t')
                #depth410660.write(str(depths[l410[maxmap410]])+'\t')
                #depth410660.write(str(depths[l660[maxmap660]])+'\t')

                
                if self.VOL.volumeweight[i, j, l410[maxmap410]] < mincoverage:
                    dmap[i, j] = 1000.

        # Prepare map
        m = Basemap(projection='merc', llcrnrlat=np.min(self.VOL.grid_lat), urcrnrlat=np.max(self.VOL.grid_lat), llcrnrlon=np.min(self.VOL.grid_lon), urcrnrlon=np.max(self.VOL.grid_lon), lat_ts=20, resolution='i')

        m.drawparallels(np.arange(0, 90, 5.), labels=[1, 0, 0, 1], linewidth=0.5, dashes=[4, 2], labelstyle='+/-', fontsize=20)
        m.drawmeridians(np.arange(-180, -110, 10.), labels=[1, 0, 0, 1], linewidth=0.5, dashes=[4, 2], labelstyle='+/-', fontsize=20)
        m.drawcountries()
        coasts = m.drawcoastlines(zorder=1, color='k', linewidth=1)

        xx, yy = np.meshgrid(self.VOL.grid_lon, self.VOL.grid_lat)
        x, y = m(xx, yy)

        #cs = plt.contourf(x, y, dmap.T, vmin=220.,levels=np.linspace(220., 280., 81.), cmap=cm.RdBu) 
        cs = plt.pcolor(x, y, dmap.T, vmin=220., vmax=280., cmap=cm.RdBu, linewidth=0, rasterized=True)


        cs.cmap.set_under([0.8, 0.8, 0.8])
        cs.cmap.set_over([0.8, 0.8, 0.8])
        cb = plt.colorbar()
        cb.set_label('Transition zone thickness (km)', size=30)
        cb.set_ticks([220,235,250,265,280])
        cb.ax.tick_params(labelsize=30)
        cb.solids.set_rasterized(True)
        xt, yt = m(-13.2, 70.6)

        m.drawcoastlines(zorder=1, color='k', linewidth=1)
        dmapall = np.ravel(dmap)

        lonmin = -158
        lonmax = -154
        latmin =  57.5
        latmax =  60
        #plot_rectangle(m, lonmin,lonmax,latmin,latmax, col="Black")

        #l = [l for l in range(len(dmapall)) if dmapall[l] > mindepth+1. and dmapall[l] < maxdepth - 1.]
        #print ('median', np.median((dmapall[l])))
        #print ('variance', np.var((dmapall[l])))




    def plot_mtzwidth_write(self,name='Megavolume',filter='rff2',conversion='EU60',factor=2.):

        #This routine is also used to make a txt file with the significant depths of the 410 and the 660
        #depth410660 = open('/raid3/annemijn/scripts/CCP/mincov40MTZ_'+conversion+'_'+filter+'_'+str(int(factor))+'.txt', 'w')
        depth410660 = open('/raid1/annemijn/scripts/CCP/MTZsign_mincov40_'+conversion+'_'+filter+'_'+str(int(factor))+'.txt', 'w')

        
        plt.figure(figsize=(18, 8))
        depths = self.VOL.grid_depth
        l410 = [x for x in range(len(depths)) if depths[x] > 380 and depths[x] < 430]
        l660 = [x for x in range(len(depths)) if depths[x] > 630 and depths[x] < 700]


        thickness1D = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        d4101D = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        d6601D = np.empty((len(self.VOL.grid_lon), len(self.VOL.grid_lat)))
        for i in range(len(self.VOL.grid_lon)):
            for j in range(len(self.VOL.grid_lat)):


                RF = self.VOL.volume[i, j,:]/self.VOL.volumeweight[i, j,:]
                std = 1.96*np.sqrt(self.VOL.volumesigma[i, j,:]/(self.VOL.volumeweight[i, j,:]*self.VOL.volumeweight[i, j,:]))
                max410 = np.argmax(RF[l410])
                max660 = np.argmax(RF[l660])
                maxamp410 = np.max(RF[l410])
                maxamp660 = np.max(RF[l660])

                # If both picks are significant, store thickness
                if RF[l410[max410]] > std[l410[max410]] and RF[l660[max660]] > std[l660[max660]]:
                    d4101D[i, j] = depths[l410[max410]]
                    d6601D[i, j] = depths[l660[max660]]
                    thickness1D[i, j] = (depths[l660[max660]]-depths[l410[max410]])

                    if self.VOL.volumeweight[i,j,l410[max410]]>=40:

                    #This routine is also used to make a txt file with the depth of the 410 and the 660
                        depth410660.write(str(self.VOL.grid_lat[j])+'\t')
                        depth410660.write(str(self.VOL.grid_lon[i])+'\t')
                        depth410660.write(str(depths[l410[max410]])+'\t')
                        depth410660.write(str(depths[l660[max660]])+'\t')
                        depth410660.write(str(maxamp410)+'\t')
                        depth410660.write(str(maxamp660)+'\n')



                    
        # Prepare map
        m = Basemap(projection='merc', llcrnrlat=np.min(self.VOL.grid_lat), urcrnrlat=np.max(self.VOL.grid_lat), llcrnrlon=np.min(self.VOL.grid_lon), urcrnrlon=np.max(self.VOL.grid_lon), lat_ts=20, resolution='i')
        m.drawparallels(np.arange(np.min(self.VOL.grid_lat), np.max(self.VOL.grid_lat), 10.), labels=[0, 0, 0, 0])#[1,0,0,1])
        m.drawmeridians(np.arange(np.min(self.VOL.grid_lon), np.max(self.VOL.grid_lon), 10.), labels=[0, 0, 0, 0])#[1,0,0,1])
        m.drawcoastlines(color='k')
        m.drawcountries(color='k')
        m.drawmapboundary(fill_color=[1.0, 1.0, 1.0])


        xx, yy = np.meshgrid(self.VOL.grid_lon, self.VOL.grid_lat)
        x, y = m(xx, yy)
        cs = plt.contourf(x, y, thickness1D.T, levels=np.linspace(220., 300., 81.), cmap=cm.gist_earth_r)
        cs.cmap.set_under('w')
        cs.cmap.set_over('w')

        plt.colorbar()
        plt.title('MTZ width')

        #This routine is also used to make a txt file with the depth of the 410 and the 660
        depth410660.close()



#
# Plot moveout aligned on 410 or 660
#

    def plot_moveout(self,d660=True,name='Megavolume',filter='rff2',conversion='EU60',factor=2.):
    # Picks all profiles in grid and organizes them by 660 depth
        plt.figure(figsize=(8, 8))
        depths = self.VOL.grid_depth

        if d660:
            ldisc = [x for x in range(len(depths)) if depths[x] > 630 and depths[x] < 730]
            disc = np.arange(630., 690., 2.)
            ldisc2 = [x for x in range(len(depths)) if depths[x] > 360 and depths[x] < 460]
        else:
            ldisc = [x for x in range(len(depths)) if depths[x] > 380 and depths[x] < 430]
            disc = np.arange(380, 440, 2.)




        weights = np.zeros((len(disc),))
        moveout = np.zeros((len(disc), len(depths)))
        d660l = []
        d410l = []
        for i in range(len(self.VOL.grid_lon)):
            for j in range(len(self.VOL.grid_lat)):
                RF = self.VOL.volume[i, j,:]
                for k in range(len(depths)):
                    if self.VOL.volumeweight[i, j, k] > 0:
                        RF[k] = RF[k]/self.VOL.volumeweight[i, j, k]
                std = 1.96*np.sqrt(self.VOL.volumesigma[i, j,:]/(self.VOL.volumeweight[i, j,:]*self.VOL.volumeweight[i, j,:]))

                maxdisc = np.argmax(RF[ldisc])


                if RF[ldisc[maxdisc]] > std[ldisc[maxdisc]] :
                    d660l.append(depths[ldisc[maxdisc]])
                    n = np.nanargmin(np.abs(depths[ldisc[maxdisc]]-disc))
                    weights[n] = weights[n]+1.
                    moveout[n,:] = moveout[n,:]+RF

                    maxdisc = np.nanargmax(RF[ldisc2])
                    d410l.append(depths[ldisc2[maxdisc]])


        d660l = np.array(d660l)
        d410l = np.array(d410l)
 
        d660c = np.arange(640., 690., 2.)
        d410c = []
        err = []
        for i in range(len(d660c)):
            filt = [d410l[x] for x in range(len(d410l)) if (d660l[x] == d660c[i] and d410l[x] > 360 and d410l[x] < 460)]
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
                    l410 = [x for x in range(len(depths)) if depths[x] > 380 and depths[x] < 440]
                    max410 = np.argmax(moveout[i, l410])
                    discp2.append(depths[l410[max410]])
                else:
                    max410 = np.argmax(moveout[i, ldisc])
                    discp.append(depths[ldisc[max410]])
                    midp.append(disc[i])
                    l660 = [x for x in range(len(depths)) if depths[x] > 630 and depths[x] < 700]
                    max660 = np.argmax(moveout[i, l660])
                    discp2.append(depths[l660[max660]])
      

        ax = plt.subplot2grid((3, 6), (0, 0), colspan=5)
 
        if d660:
            n, bins, patches = plt.hist(d660l, np.arange(639., 713., 2.), histtype='bar')
            data = zip(bins+1., n) # +1 as the stack resolution is only at even values
            cmbl = plt.cm.get_cmap('bone_r')
            for d, p in zip(n, patches):
                plt.setp(p, 'facecolor', cmbl(d/max(n)*0.6+0.2))
            # f=open('histogram_660.txt','wb')
            # np.savetxt(f,data,delimiter='\t')
        else:
            plt.hist(d660l, np.arange(389, 431, 2.), histtype='bar')
        plt.ylabel('# of grid points')
        plt.text(645, 600, 'a.', fontsize=16)
        if d660:
            # plt.xlabel('660 depth (km)')
            plt.xlim([644.9, 691.07])
        else:
            plt.xlabel('410 depth (km)')
            plt.xlim([390, 430])
        plt.gca().set_yticks([0, 100, 300, 500, 700])
        plt.gca().set_xticks([650, 660, 670, 680, 690])
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
            plt.xlim([644.9, 691.07])
            plt.xlabel('660 depth (km)')
        else:
            plt.xlim([390, 430])
            plt.xlabel('410 depth (km)')
        plt.ylim([350, 750])
        plt.gca().invert_yaxis()
        plt.text(646, 380, 'b.', fontsize=16)
        box = plt.gca().get_position()
        axcb = plt.axes([box.x0*1.05+box.width*1.05, box.y0+0.1, 0.01, box.height-0.2])
        cb = plt.colorbar(cs, cax=axcb, orientation='vertical')
        cb.set_ticks([-0.1, 0., .1])
        cb.set_label('relative amplitude')

        plt.show()



def plot_rectangle(bmap, lonmin,lonmax,latmin,latmax, col):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    bmap.plot(xs, ys,latlon = True, color=col)

def plot_line (bmap,lonmin,lonmax,latmin,latmax,col):
    xs = [lonmax,lonmax,lonmin]
    ys = [latmin,latmin,latmax]
    bmap.plot(xs, ys,latlon = True, color="Black")
