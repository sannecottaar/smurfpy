# -*- coding: iso-8859-1 -*-
'''
implementation of AFR_RF_CR1 combined Receiver function and crustal model for RF depth correction
Model is detailed in supplement of African tomography model AFRP20 outlined in Boyce et al., 2020 submitted to Gcubed
The third velocity value is the vp from the 50km depth interval within AFRP20, converted to vs as well.
'''

import numpy as np


import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap

import math
import scipy

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

class AFR_RF_CR1_MOD_model(object):
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
    
    self.depthbounds=[]
    self.depth=[]
    self.Vp=[]
    self.Vs=[]
    self.topo=[]
#    self.dVp=[]
#    self.dVs=[]



  #########################################################################
  #- read a 3D model
  #########################################################################

  def read(self,directory='./', verbose=False):
    """ read an ses3d model from a file

    read(self,directory,filename,verbose=False):
    """

    #- read block files ====================================================

    lat=[]
    lon=[]
    depthbounds=[]
    vp =[]
    vs=[]
    topo=[]
#    This file is taken from matlab crustal generation script.
#        Topo currently taken to be zero
#        Lat Lon TOPO D1 D2 Vp1 Vp2 AFRP20_vp_50km Vs1 Vs2  AFRP20_vs_50km
    fid=open(directory+'AFR_RF_deps_ak135_vels_AFRP20_RF_CR1_mod.txt','r')

    for line in fid.readlines():
        val=line.split()
        lat.append(float(val[0]))
        lon.append(float(val[1]))
        topo.append(float(val[2])) # Make sure this topography is going is as negative. Depths are positive.
        depthbounds.append([float(val[2]), float(val[3]), float(val[4])])
        vp.append([float(val[5]),float(val[6]),float(val[7])])
        vs.append([float(val[8]),float(val[9]),float(val[10])])
    
    depthboundst=np.array(depthbounds)
    depthbounds=depthboundst.astype(np.float)
    vpt=np.array(vp)
    vp = vpt.astype(np.float)
    vst=np.array(vs)
    vs=vst.astype(np.float)
    topot = np.array(topo)
    topo=topot.astype(np.float)
    
    self.lon=np.unique(lon)
    self.lat=np.unique(lat)
    
    self.depthbounds=depthbounds.reshape((len(self.lat),len(self.lon),3),order='C')
    self.Vp=vp.reshape((len(self.lat),len(self.lon),3),order='C')
    self.Vs=vs.reshape((len(self.lat),len(self.lon),3),order='C')
    self.topo = topo.reshape((len(self.lat),len(self.lon)),order='C')

    self.lon_min=np.min(self.lon)
    self.lon_max=np.max(self.lon)
    self.lat_min=np.min(self.lat)
    self.lat_max=np.max(self.lat)
    self.dlon=self.lon[1]-self.lon[0]
    self.dlat=self.lat[1]-self.lat[0]

  #########################################################################
  #- plot horizontal slices
  #########################################################################

  def plot_slice(self,depth,whattoplot='Vp',colormap='tomo',res='i',latlim=[-90,90],lonlim=[-180,180],verbose=False):
    """ plot horizontal slices through an ses3d model

    plot_slice(self,depth,colormap='tomo',res='i',verbose=False)

    depth=depth in km of the slice
    colormap='tomo','mono'
    res=resolution of the map, admissible values are: c, l, i, h f

    """
 
    #- set up a map and colourmap -----------------------------------------
    if latlim and lonlim:
        m=Basemap(projection='merc',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1],lat_ts=20,resolution=res)
    else:
        m=Basemap(projection='merc',llcrnrlat=self.lat_min,urcrnrlat=self.lat_max,llcrnrlon=self.lon_min,urcrnrlon=self.lon_max,lat_ts=20,resolution=res)
    m.drawparallels(np.arange(self.lat_min,self.lat_max,10.),labels=[1,0,0,1])
    m.drawmeridians(np.arange(self.lon_min,self.lon_max,10.),labels=[1,0,0,1])

    m.drawcoastlines()
    m.drawcountries()

    m.drawmapboundary(fill_color=[1.0,1.0,1.0])

    lon=np.arange(lonlim[0],lonlim[1],.5)
    lat=np.arange(latlim[0],latlim[1],.5)
    

    xx,yy= np.meshgrid(lon,lat)

    x,y=m(xx.T,yy.T)
    toplot=np.zeros_like(x)
    for i in range(len(lon)):
        for j in range(len(lat)):
            toplot[i,j]=self.get_value(depth,lon[i],lat[j],whattoplot)
            

#    print(toplot)
    minval=np.min(np.min(toplot))
    maxval=np.max(np.max(toplot))
#    print(minval,maxval)
#    ref=np.mean(np.mean(toplot))
    contours=np.round(np.linspace(minval,maxval,10),2)

    plt.contourf(x,y,toplot,contours,cmap=plt.cm.get_cmap('RdBu'))
    plt.colorbar()
    plt.title(whattoplot+' at '+str(depth)+' km')

    ################################
    # get value at specific location
    ################################
  def get_value(self,depth,lon,lat,whatmodel='Vp',method='nearest'):

      if method=='nearest': # find nearest neighbor in lat and lon. Interpolate in depth
          if lon>self.lon_min-0.51 and lon<self.lon_max+0.51 and lat>self.lat_min-0.51 and lat < self.lat_max+0.51:
              #- loop over subvolumes to collect information ------------------------

                
              lonind=np.argmin(np.abs(self.lon-lon))
              latind=np.argmin(np.abs(self.lat-lat))
              
              if whatmodel == 'topo':

                  return self.topo[latind,lonind]
              else:
                  depths=self.depthbounds[latind,lonind,:]
    
                  ind = [x for x in range(len(depths)) if depths[x]<=depth][-1]
    
                  model=getattr(self,whatmodel)
                  val=model[latind,lonind,ind]
    
                  return val
          else:
              print('value requested outside for ', lat, lon, depth)




if __name__ == "__main__":

    
    modcrust=AFR_RF_CR1_MOD_model()
    modcrust.read(directory='./', verbose=True)
    print(modcrust.get_value(50,30.23,-20.44,whatmodel='topo'))
    plt.figure(figsize=(6,6))
    # Note: the mask for this model is not implemented
    modcrust.plot_slice(50.,whattoplot='Vp',latlim=[-40,40],lonlim=[-20,60])
    plt.show()  
