# -*- coding: iso-8859-1 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 10
from mpl_toolkits.basemap import Basemap
from geographiclib.geodesic import Geodesic as geo

import math
import scipy
import scipy.ndimage
from scipy import interpolate


#########################################################################
#- define ses3d model class
#########################################################################

def haversine(lat1, long1, lats2, longs2):
    """
    Calculate the distance between two points on earth in m
    """
    d=[]

    for i in range(len(lats2)):
        lat2=lats2[i]
        long2=longs2[i]
        earth_radius = 6371.e3  # m
        dLat = math.radians(lat2 - lat1)
        dLong = math.radians(long2 - long1)

        a = (math.sin(dLat / 2) ** 2 +
             math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dLong / 2) ** 2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        d.append(earth_radius * c)
    return float(d[0])


def godown(v1,rad1,v2,rad2,slow):
    # calculates time (s) and radial distance (rad) for one segment
    dr=rad1-rad2
    v_av=(v1+v2)/2.
    slow=slow*(180./(np.pi))#!! Slowness is given in sec/deg -- convert to sec/rad
    dtheta=slow*dr/rad1*1/np.sqrt((rad1/v_av)**2-slow**2)
    dtime=dr/rad1*(rad1/v_av)**2/np.sqrt((rad1/v_av)**2-slow**2)
    return dtheta,dtime

class AFRP20_RF_CR1_model(object):
    """ class for reading, writing, plotting and manipulating and ses3d model
    """
    
    def __init__(self):
        """ initiate the ses3d_model class
        initiate list of submodels and read rotation_parameters.txt
        """
    
        self.nsubvol=0
        self.lat_min=0.0
        self.lat_max=0.0
        self.lon_min=0.0
        self.lon_max=0.0
        self.lat_centre=0.0
        self.lon_centre=0.0
    
        self.lat=[]
        self.lon=[]
    
        self.depth=[]
    #    self.Vp=[]
    #    self.Vs=[]
        self.dVp=[]
        self.dVs=[]
    
    
    
      #########################################################################
      #- read a 3D model
      #########################################################################

    def read(self,directory,filenames=None,verbose=False):
        """ read an ses3d model from a file
    
        read(self,directory,filename,verbose=False):
        """
    
        #- read block files ====================================================
        lon=[]
        lat=[]
        depth=[]
        vp=[]
        vs=[]
    
        fid=open(directory+filenames[0],'r')
        for line in fid.readlines():
            val=line.split()
            lon.append(float(val[2]))
            lat.append(float(val[1]))
            depth.append(float(val[0]))
            vp.append(float(val[3]))
    
    #    sorted W->E, S->N, C->CMB
        lon=np.array(lon) #-360 # use (-360) to correct to -180/180, already in correct format.
        lat=np.array(lat)
        depth=np.array(depth)
        vp = np.array(vp)
        vs = vp*(depth/2891.+2.) # from Ritsema et al., (2011)
        
        self.lon=np.unique(lon)
        self.lat=np.unique(lat)
        self.depth=np.unique(depth)
    
        self.dVp=vp.reshape((len(self.depth),len(self.lat),len(self.lon)),order='C')
        self.dVs=vs.reshape((len(self.depth),len(self.lat),len(self.lon)),order='C')
        # flipping latitudes
        # swaps place of 2nd and 3 dimension.
        self.dVp=np.transpose(self.dVp,(0,2,1))
        self.dVs=np.transpose(self.dVs,(0,2,1))
    
    #    self.RVpVs=self.dVp/self.dVs
    
        self.lon_min=np.min(self.lon)
        self.lon_max=np.max(self.lon)
        self.lat_min=np.min(self.lat)
        self.lat_max=np.max(self.lat)
        self.dlon=self.lon[1]-self.lon[0]
        self.dlat=self.lat[1]-self.lat[0]

  #########################################################################
  #- plot horizontal slices
  #########################################################################

    def plot_slice(self,depth,whattoplot='dVp',colormap='tomo',res='h',latlim=None,lonlim=None,verbose=False):
        """ plot horizontal slices through an ses3d model
    
        plot_slice(self,depth,colormap='tomo',res='i',verbose=False)
    
        depth=depth in km of the slice
        colormap='tomo','mono'
        res=resolution of the map, admissible values are: c, l, i, h f
    
        """
     
        #- set up a map and colourmap -----------------------------------------
        if latlim and lonlim:
            m=Basemap(projection='merc',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1],lat_ts=0,resolution=res)
        else:
    #        m=Basemap(projection='merc',llcrnrlat=self.lat_min,urcrnrlat=self.lat_max,llcrnrlon=self.lon_min,urcrnrlon=self.lon_max,lat_ts=0,resolution=res)
            m=Basemap(projection='merc',llcrnrlat=-40,urcrnrlat=40,llcrnrlon=-20,urcrnrlon=60,lat_ts=0,resolution=res)
    #        m=Basemap(projection='merc',llcrnrlat=5,urcrnrlat=15,llcrnrlon=38,urcrnrlon=48,lat_ts=0,resolution=res)
    
        
        m.drawparallels(np.arange(self.lat_min,self.lat_max,10.),labels=[1,0,0,1])
        m.drawmeridians(np.arange(self.lon_min,self.lon_max,10.),labels=[1,0,0,1])
    
        m.drawcoastlines()
        m.drawcountries()
    
        m.drawmapboundary(fill_color=[1.0,1.0,1.0])
    
        model=getattr(self,whattoplot)
        layer=np.argmin(np.abs(self.depth-depth))
        print('True plotting depth is ', self.depth[layer])
    
        xx,yy= np.meshgrid(self.lon,self.lat)
        x,y=m(xx,yy)
        x=x.T
        y=y.T
     
        minval=np.min(np.min(model[layer,:,:]))
        maxval=np.max(np.max(model[layer,:,:]))
        print(maxval)
        print(minval)
        if minval<-6:
            minval=-6
        if maxval>6:
            maxval=6
    #    ref=np.mean(np.mean(model[layer,:,:]))
        contours=np.round(np.linspace(-2,2,9),2)
        toplot=model[layer,:,:]
        
    #    print(x.shape)
    #    print(y.shape)
    #    print(toplot.shape)
        
        
        plt.contourf(x,y,toplot,contours,cmap=plt.cm.get_cmap('RdBu'), extend='both')
        plt.colorbar()
    
        plt.cm.get_cmap('RdBu').set_over('blue')
        plt.cm.get_cmap('RdBu').set_under('red')
        
        plt.title(whattoplot+' at '+str(depth)+' km')


    def plot_XC(self,lon1=18,lon2=48,lat1=-38,lat2=25,numpoints=200,Upper_Depth=0, Lower_Depth=700,whattoplot='dVp', P_lim=2):
    
        ''' Function to plot cross sections through the tomographic model. '''
        
        # Sort out the coross section parameters
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

        dist = []
        for i in range(len(lats)):
            dist.append(haversine(lats[0], lons[0], [lats[i]], [lons[i]])/111194.)    

        # Now deal with the model volume itself.
        row = (lons-np.min(self.lon))/(self.lon[1]-self.lon[0])
        for i in range(len(row)):
            if row[i] < 0:
                row[i] = row[i]+len(self.lon)
        col = (lats-np.min(self.lat))/(self.lat[1]-self.lat[0])
        
        model=getattr(self,whattoplot)
        crossec = []
        for dp in range(len(self.depth)):
            #print('True plotting depth is ', self.depth[dp])
            #print(dp)
            #print(np.shape(model[:,:,dp]))
            crossec.append(scipy.ndimage.map_coordinates(model[dp,:,:], np.vstack((row, col))))    
        
        crossec = np.array(crossec)
        depths = self.depth
        xx, yy = np.meshgrid(dist, depths)
        #print(depths)
        #print(dist)
        
        #print(np.shape(crossec))
        #print(np.max(np.max(crossec)))
        #print(np.min(np.min(crossec)))


#        Interpolate the data to facilitate smooth plotting.
        f = interpolate.interp2d(dist, depths, crossec, kind='quintic')
        
    
        dist_int=0.1
        depth_int=10
        dist_interp=np.arange(np.min(dist), np.round(np.max(dist),2), dist_int)
        depths_interp=np.arange(np.min(depths), np.round(np.max(depths),2), depth_int)
        
        xx_i, yy_i = np.meshgrid(dist_interp, depths_interp)

        crossec_i=f(dist_interp,depths_interp)
        
        plot_width=7
        plot_height=2.5
        
        fig = plt.figure(figsize =(plot_width,plot_height))
        
        ax1 = fig.add_axes([0.1, 0.2, 0.85, 0.7], projection=None, polar=False,facecolor='white',frame_on=True)
        ax1.set_ylabel('Depth (km)', fontsize=12)
        ax1.set_xlabel('Angular distance (dg)', fontsize=12)
        ax1.set_ylim([Upper_Depth, Lower_Depth])
        ax1.invert_yaxis()
        ax1.set_xticks(np.arange(0,max(dist),10))
        ax1.set_yticks(np.arange(Upper_Depth, Lower_Depth+1,200))
        ax1.tick_params(labelsize=12)
        
        #Plot axes for gaps between XCs
        axgb = fig.add_axes([0.1, 0.0, 0.85, 0.1], projection=None, polar=False,facecolor='white',frame_on=False)
        axgb.set_xticks([])
        axgb.set_yticks([])

        #Plot axes for gaps between XCs
        axgt = fig.add_axes([0.1, 0.9, 0.85, 0.1], projection=None, polar=False,facecolor='white',frame_on=False)
        axgt.set_xticks([])
        axgt.set_yticks([])
    
        
        contours=np.round(np.linspace(-P_lim,P_lim,40+1),2)
        #im1=ax1.contourf(xx,yy,crossec,contours,cmap=plt.cm.get_cmap('RdBu'), extend='both')
        im1=ax1.contourf(xx_i,yy_i,crossec_i,contours,cmap=plt.cm.get_cmap('RdBu'), extend='both')

        #im1=ax1.pcolormesh(xx, yy, crossec, vmin=-4, vmax=4, rasterized=False,cmap=plt.cm.get_cmap('RdBu'))
        fig.colorbar(im1,ax=ax1, shrink=0.8)
        
        left_label = axgt.text(0, 0.8, str(lon1) + 'E,'+str(lat1)+'N', ha="center",va="top",fontsize=10, color='black', visible=True, clip_on=False, bbox=dict(facecolor='white',edgecolor='black', pad=2.0)) 

        right_label = axgt.text(0.8, 0.8, str(lon2) + 'E,'+str(lat2)+'N', ha="center",va="top",fontsize=10, color='black', visible=True, clip_on=False, bbox=dict(facecolor='white',edgecolor='black', pad=2.0)) 

        
    
 #########################################################################
  #- retrieve horizontal slices
  #########################################################################

    def get_slice(self,depth,whattoplot='dVp',colormap='tomo',res='i',verbose=False):
        """ plot horizontal slices through an ses3d model
    
        plot_slice(self,depth,colormap='tomo',res='i',verbose=False)
    
        depth=depth in km of the slice
        colormap='tomo','mono'
        res=resolution of the map, admissible values are: c, l, i, h f
    
        """
    
        model=getattr(self,whattoplot)
        layer=np.argmin(np.abs(self.depth-depth))
        xx,yy= np.meshgrid(self.lon,self.lat)
     
    
        return xx.T, yy.T, model[layer,::,:]  

    ################################
    # get value at specific location
    ################################
    def get_value(self,depth,lon,lat,whatmodel='dVp',method='nearest'):
 
        if method=='nearest': # find nearest neighbor in lat and lon. Interpolate in depth
            if lon>self.lon_min and lon<self.lon_max and lat>self.lat_min and lat < self.lat_max and depth <= np.max(self.depth):
                #- loop over subvolumes to collect information ------------------------
                lonind=np.argmin(np.abs(self.lon-lon))
                latind=np.argmin(np.abs(self.lat-lat))
                model=getattr(self,whatmodel)
                valloc=model[:,lonind,latind]
                val = np.interp(depth,self.depth,valloc)
                return val
            else:
                print('value requested outside for ', lat, lon, depth)
                return 0
            


if __name__ == "__main__":
    
    mod=AFRP20_RF_CR1_model()
    mod.read('../MODELS/AFRP20_RF_CR1/', filenames = ['AFRP20_RF_CR1_P_model.dat'], verbose=True)

#    Plotting
    #plt.figure(figsize=(6,6))
    #mod.plot_slice(100.,whattoplot='dVp')
    
    UD=50; LD=800
    mod.plot_XC(lon1=21,lon2=44,lat1=-34,lat2=17,numpoints=120,Upper_Depth=UD,Lower_Depth=LD,whattoplot='dVp',P_lim=4)
    plt.savefig('./FIGURES/AFRP20CR_XC_dVp_'+str(UD)+'-'+str(LD)+'km.pdf')

##    Test get model value
#    dVp=mod.get_value(410,38,8,'dVp')
#    print('dVp 410: ' + str(dVp))
#
#    dVp=mod.get_value(660,38,8,'dVp')
#    print('dVp 660: ' + str(dVp))
#    
#    
#    dVs=mod.get_value(410,38,8,'dVs')
#    print('dVs 410: ' + str(dVs))
#
#    dVs=mod.get_value(660,38,8,'dVs')
#    print('dVs 660: ' + str(dVs))  
#    
    
    
#    Median velocity at each depth interval
#    for d in mod.depth:
#        x,y, slicet = mod.get_slice(d,whattoplot='dVp')
#        #print(slicet)
#        print(d, np.nanmedian(slicet))

 
    #plt.show()
