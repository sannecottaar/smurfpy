# common conversion point stacking

# import modules
import sys
#sys.path.insert(1,'/raid2/sc845/Python/basemap-1.0.7/lib/')
#sys.path.append('/raid2/sc845/Python/basemap-1.0.7/lib/python')

#sys.path.append('/raid2/sc845/Python/lib/python/')
sys.path.append('/raid2/sc845/Tomographic_models/EU60/')
sys.path.append('/raid2/sc845/Python/geographiclib-1.34/')
from geographiclib.geodesic import Geodesic as geo
import numpy as np
#from matplotlib.mlab import griddata
import scipy
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

import mpl_toolkits





# definition of the half width of the fresnel zone
knotspacing = lambda r, vs: 1./2.*np.sqrt(((10./3.*vs)+r)**2.-r**2.) # in m for a 10s wave


def haversine(lat1, long1, lats2, longs2,depth):
    """
    Calculate the distance between two points in earth in m
    """
    d=[]

    for i in range(len(lats2)):
        lat2=lats2[i]
        long2=longs2[i]
        radius = 6371.e3 - depth*1.e3  # m
        dLat = math.radians(lat2 - lat1)
        dLong = math.radians(long2 - long1)

        a = (math.sin(dLat / 2) ** 2 +
             math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dLong / 2) ** 2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        d.append(radius * c)
    return float(d[0])


def weight(distance,depth,vs,factor):
    '''
    Calculates the weight based on the distance and the fresnel zone width for a given depth
    '''
    #print knotspacing(depth*1.e3,vs)
    delta=distance/(factor*knotspacing(depth*1.e3,vs)) # distance in m ~~~~~~ fresnel zone times factor~~~ 
    #print delta
    if delta>2:
        weight=0
    elif delta>1:
        weight=.25*(2.-delta)**3.
    else:
        weight=.75*delta**3.-1.5*delta**2.+1.
 
    return weight

class VOL(dict):
    '''
    Initiates a dictionary for the CCP volume
    '''
    def __init__(self,*arg,**kw):
        super(VOL,self).__init__(*arg,**kw)
        self.__dict__=self
    def __getattr__(self,name):
        return self[name]

class ccp_volume(object):
    """
       Handling large stacked volumes
    """
    def __init__(self,*arg,**kw):
        self.VOL=VOL(*arg,**kw)




#############################################################################
# Start a new volume dictionary
#############################################################################

    def start_empty_volume(self,name='Megavolume',filter='rff2',conversion='prem',  factor=1., lonmin=None, lonmax=None, lonrez=None, latmin=None,latmax=None,latrez=None, depmin=None,depmax=None,deprez=None):
        '''
        Start the empty CCP volume

        Parameters
        ---------
        name    String that defines the name of the directory where the volumes is stored
        filter  Sring that refers to the name of the receiver function to be used (If you add a new one, this needs to be implemented in the cod)
        conversion String for the  conversion model to be used (default = prem)
                Receiver functions needs to contain dictionary with this conversion

        factor  This sets the amount of smoothing beyond the Fresnel Zone. Recommended to start with 2 (and maybe try 1 later)
        lonmin, lonmax, lonrez Longitude bounds and number of steps to set up the grid
        latmin, latmax, latrez Latitude bounds and number of steps to set up the grid
        depmin, depmax, deprez Depths bounds and number of steps to set up the grid
        '''

        #Make directory for the volumes
        dirout='../CCP_volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)
        if not os.path.exists(dirout):
            os.makedirs(dirout)
        if not os.path.exists(dirout+'/RF_lists/'):
            os.makedirs(dirout+'/RF_lists/')
        # Store initial empty PICKLE
        outfilename='../CCP_volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/Stack_0.PICKLE'
        donefilename='../CCP_volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/list_of_stacks.txt'

        #setup volume grid
        grid_depth=np.linspace(depmin,depmax,num=deprez)
        grid_lon=np.linspace(lonmin,lonmax,num=lonrez)
        grid_lat=np.linspace(latmin,latmax,num=latrez)


        # Read in STW105 for the velocities used in the fresnel zone width
        # (I forgett why I use this model instead of prem. Not that it matters much...it's only used for the weighting factors)
        #Replace with PREM
        table=[] 
        for line in open('../Tools/STW105.txt').readlines()[3:]:
            if line[0]!='#':
                numbers= list(map(float,line.split()))
                table.append(numbers)
        table=np.array(table)
        refdepths=(6371.e3-table[:,0])
        refVp=(table[:,6])
        refVs=(table[:,7])

        ####### Define volumes
        volume=np.zeros([len(grid_lon),len(grid_lat),len(grid_depth)])
        volumeweight=np.zeros([len(grid_lon),len(grid_lat),len(grid_depth)]) #tracks summed weight
        volumesigma=np.zeros([len(grid_lon),len(grid_lat),len(grid_depth)]) # tracks weighted difference
        volumesign=np.zeros([len(grid_lon),len(grid_lat),len(grid_depth)]) #tracks stack of sign
        num=np.zeros([len(grid_lon),len(grid_lat),len(grid_depth)]) # tracks number of RFs involved

        grid_vs= np.interp(grid_depth,refdepths[::-1]/1.e3,refVs[::-1])
        grid_vp= np.interp(grid_depth,refdepths[::-1]/1.e3,refVp[::-1])

        # Store values and write out volume
        self.VOL['lonmin']=lonmin
        self.VOL['lonmax']=lonmax
        self.VOL['latmin']=latmin
        self.VOL['latmax']=latmax
        self.VOL['depmin']=depmin
        self.VOL['depmax']=depmax
        self.VOL['grid_depth']=grid_depth
        self.VOL['grid_lon']=grid_lon
        self.VOL['grid_lat']=grid_lat
        self.VOL['volume']=volume
        self.VOL['volumeweight']=volumeweight
        self.VOL['volumesign']=volumesign
        self.VOL['volumesigma']=volumesigma
        self.VOL['num']=num
        self.VOL['count']=0
        self.VOL['grid_vs']=grid_vs
        self.VOL['grid_vp']=grid_vp



        with open(outfilename,'wb') as handle:
            msgpack.pack(self.VOL,handle)
            handle.close()

        with open('../CCP_volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/filenames.dat','w') as handle:
            handle.write('%d %s \n'% (0, outfilename))
            handle.close()


        print('DONE CREATING EMPTY VOLUME')
        return self




#############################################################################
# Load latest volume to dictionary
#############################################################################
    def load_latest(self,name='Megavolume',filter='rff2',conversion='prem',factor=1.):
        '''
        Loads latest volume
        '''
        print(name)
        line=open('Volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/filenames.dat','r').readlines()[-1]
        runnum=int(float(line.split()[0]))
        volumefile=line.split()[1]
        print(runnum, volumefile)



        # get last stackfile name

        ####### Read in volumes

        self.VOL.update(msgpack.unpack(open(volumefile,'rb'), use_list=False,object_hook=m.decode))
        #del self.VOL['trackRFs']

        
        return self


#############################################################################
# Add list of receiver functions to Volume
#############################################################################

    def addlist(self,rflist,name='Megavolume',filter='rff2',conversion='EU60',factor=1.):
        ''' 
        Adds list of PICKLE files to volume
        The code has no trouble adding the same PICKLE file twice, which means those are not double weighted in the stack. Might need a fix...
        '''
        ## set volume lats and lons
        rffilter=filter
        line=open('../CCP_volumes/'+name+'_'+rffilter+'_'+conversion+'_'+str(factor)+'/filenames.dat','r').readlines()[-1]
        runnum=int(float(line.split()[0]))
        volumefile=line.split()[1]
        print(runnum, volumefile)
 
        rffile=open('../CCP_volumes/'+name+'_'+rffilter+'_'+conversion+'_'+str(factor)+'/RF_lists/rflist'+(str(runnum+1))+'.dat','a')

        for i in range(len(rflist)): #range(cat.count()): # loop through all receiver functions for a given station
            print(i, 'out of', len(rflist))
            if os.path.isfile(rflist[i]):
                #try:
                    # Read and retrieve receiver function
                    seis=read(rflist[i],format='PICKLE')
                    self.VOL.count=self.VOL.count+1
                    rffile.write("%d %s \n" % (int(self.VOL.count),rflist[i]))


                    RF = np.real(getattr(seis[0],filter)['iterativedeconvolution'])
                  
                    indm=np.argmax(np.abs(RF))   
                    if np.mean(RF[indm-10:indm+10])<0.: # flip receiver function if needed
                        RF=RF*-1.
                    RF=RF/np.max(np.abs(RF)) # Normalize RF
                    for d in range(len(self.VOL.grid_depth)): # loop through all depths
                            # find (lat,lon) of the ray path at this given depth 3D
                            x       = np.argmin(np.abs((seis[0].conversions[conversion]['depths']-self.VOL.grid_depth[d])))
                            #find (lat,lon) of the ray path at this given depth 3D
                            latx = seis[0].conversions[conversion]['latitudes'][x]
                            lonx = seis[0].conversions[conversion]['longitudes'][x]
                            x_RF    = np.argmin(np.abs((seis[0].conversions[conversion]['depthsfortime']-self.VOL.grid_depth[d]))) 

                            # loop through all lats and lons for the given depth
                            # put some limits on grid to use
                            lonind  = np.argmin(np.abs(self.VOL.grid_lon-lonx))
                            inds=int(round(d/30.)+3.) # widen as going deeper
                            lonlim  = np.arange(lonind-inds,lonind+inds)
                            latind  = np.argmin(np.abs(self.VOL.grid_lat-latx))
                            latlim  = np.arange(latind-inds,latind+inds)
                            for k in lonlim:
                                for j in latlim:
                                  if k>-1 and j>-1 and k < len(self.VOL.grid_lon) and j < len(self.VOL.grid_lat):
                                    # Stack into 1D stacks
                                    dist = haversine(latx,lonx,[self.VOL.grid_lat[j]],[self.VOL.grid_lon[k]],self.VOL.grid_depth[d]) # calculate distance

                                    w    = weight(dist,self.VOL.grid_depth[d],self.VOL.grid_vs[d],factor) # calculate weight using S wave fresnel zone
                                    if w>0:
                                        #self.VOL.trackRFs[k][j][d].append(self.VOL.count)
                                        self.VOL.num[k,j,d]          = self.VOL.num[k,j,d]+1.                             # count number of receiver functions

                                        self.VOL.volume[k,j,d]       = self.VOL.volume[k,j,d]+w*RF[x_RF]                  # stack receiver function into volume  
                                        self.VOL.volumeweight[k,j,d]    = self.VOL.volumeweight[k,j,d]+w                     # stack weights    
                                        self.VOL.volumesigma[k,j,d]  = self.VOL.volumesigma[k,j,d]+w*(RF[x_RF]-(self.VOL.volume[k,j,d]/self.VOL.volumeweight[k,j,d]))**2. # stack sign 
                                        self.VOL.volumesign[k,j,d]   = self.VOL.volumesign[k,j,d]+w*np.sign(RF[x_RF])     # stack sign of receiver function

        print('RFs used', self.VOL.count)

        # Write out at the end of list
        outfilename='../CCP_volumes/'+name+'_'+rffilter+'_'+conversion+'_'+str(factor)+'/Stack_'+str(int(runnum+1))+'.PICKLE'
        with open(outfilename,'wb') as handle:
                    msgpack.pack(self.VOL,handle)
                    handle.close()
        with open('../CCP_volumes/'+name+'_'+filter+'_'+conversion+'_'+str(factor)+'/filenames.dat','a') as handle:
            handle.write('%d %s \n'% (runnum+1, outfilename))
            handle.close()
        rffile.close()

        print('DONE')
        return self


