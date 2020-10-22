#
import sys

sys.path.append('../Tools/PLOTTING/')
import Africa_AFRP20_RF_CR1

sys.path.append('../Tools/MODELS/AFR_RF_CR1_MOD/')
import AFR_RF_CR1_MOD

from obspy import read
from obspy.taup import TauPyModel
import time
import glob

import numpy as np
from geographiclib.geodesic import Geodesic as geo
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

ak135=True
rffilter='jgf1'

if ak135:
    taupmodel = TauPyModel(model='ak135') # Could change to AK135 to setup in accordance with AFRP20


def godown(v1,rad1,v2,rad2,slow):
    # formulas can be found in Lay and Wallace P92.
    # calculates time (s) and radial distance (rad) for one segment
    dr=rad1-rad2
    v_av=(v1+v2)/2.
    slow=slow*(180./(np.pi))#!! Slowness is given in sec/deg -- convert to sec/rad
    dtheta=slow*dr/rad1*1/np.sqrt((rad1/v_av)**2-slow**2)
    dtime=dr/rad1*(rad1/v_av)**2/np.sqrt((rad1/v_av)**2-slow**2)
    return dtheta,dtime

def init_model():
    # Read crustal and mantle models
    global mod, crust
    mod=Africa_AFRP20_RF_CR1.AFRP20_RF_CR1_model()
    mod.read(directory='../Tools//MODELS/AFRP20_RF_CR1/', filenames = ['AFRP20_RF_CR1_P_model.dat'], verbose=True)
    crust=AFR_RF_CR1_MOD.AFR_RF_CR1_MOD_model()
    crust.read(directory='../Tools//MODELS/AFR_RF_CR1_MOD/', verbose=True)
    
def compute_onestation(dir, filter='jgf1'):

    # loop through receiver functions
#    stalist=sorted(glob.glob(dir+'*.*/*.PICKLE'))
    
    station_file = open(dir + '/selected_RFs_'+str(filter)+'.dat','r')
    stalist=station_file.read().strip().split()
    
    if len(stalist)>0:
        for i in range(len(stalist)): 
            print(i, stalist[i])
            seis=read(stalist[i],format='PICKLE') # read in receiver function

            if not 'AFRP20CR' in seis[0].conversions:
#                dt=seis[0].stats['delta']                            # set sampling read
                dist=seis[0].stats['dist']                     # set distance
                depth= seis[0].stats['evdp'] # set depth
                slat=float(seis[0].stats['stla'])
                slon=float(seis[0].stats['stlo'])
                stel=round(float(seis[0].stats['stel']),3)
#                azi=float(seis[0].stats['az'])
                bazi=float(seis[0].stats['baz'])

#                elat= seis[0].stats['evla']
#                elon= seis[0].stats['evlo']
#                times=seis[0].jgf1['time']                            # set time axis
                line=geo.WGS84.Line(slat,slon,bazi)
                label='P'
                # get ak135 reference for P wave
                arr = taupmodel.get_travel_times(depth,dist,[label])
#                Pt= arr[0].time
                Pslow=arr[0].ray_param*np.pi/180.
#                Pi=arr[0].incident_angle


                start=time.time()

                # initialize values
                xs2=[0.]        # distance Pds phase
                xp2=[0.]        # distance P phase
                tsum2=[0.]      # total time difference Pds and P (including horizontal distance)
                tpsum2=[0.]     # integrated time P phase
                tssum2=[0.]     # integrated time Pds phase
                dtcorr2=[0.]    # horizontal distance to time
                latP=[slat]     # latitude P phase
                lonP=[slon]     # longitude P phase
                latPds=[slat]   # latitude Pds phase
                lonPds=[slon]   # longitude Pds phase


                # Starting depth from station elevation and make it negative -> depth
                starting_depth=-1.0*stel
		
				# This give the topography starting at an integer value.
                #depthrange=np.arange(int(starting_depth),1500,2)

                if starting_depth < 0.0:
                    depthrange=np.arange(0.0,1500.0,2.0)
                    depthrange=np.insert(depthrange, 0, starting_depth)
                elif starting_depth == 0.0:
                    depthrange=np.arange(0.0,1500.0,2.0)
                elif starting_depth > 0.0:
                    array_start = 2*((starting_depth // 2)+1)
                    depthrange=np.arange(float(array_start),1500.0,2.0) 
                    depthrange=np.insert(depthrange, 0, starting_depth)
                else:
                    print('ERROR - something has gone quite wrong here....')
                    print('Exiting....')
                    sys.exit()

                
                print('The starting depth is: '+str(depthrange[0]))

                modrad2=6371.-depthrange
                # get slownesses for converted phases 
                Pdslow=[Pslow]
                f= interp1d(seis[0].conversions['ak135']['depths'],seis[0].conversions['ak135']['slowness']+Pslow, fill_value='extrapolate')
                for d in range(1,len(depthrange)):
                    Pdslow.append(f(depthrange[d]))
    
    
                
                topo=crust.get_value(0, slon,slat,whatmodel='topo') # This will be a negative depth
                # startin
                modvs2=[crust.get_value(topo,lonP[-1],latP[-1],'Vs')]
                modvp2=[crust.get_value(topo,lonP[-1],latP[-1],'Vp')]
 
                print('starting vel are', modvp2, modvs2, 'at', starting_depth)
                
                for d in range(1,len(depthrange)):
                    #!!! refind vp and vs
                    if depthrange[d] < topo:
                        # Station is actually at a greater elevation than that given by interpolated elevation in CRUST1.0
                        modvs2.append(crust.get_value(topo,lonP[-1],latP[-1],'Vs'))
                        modvp2.append(crust.get_value(topo,lonP[-1],latP[-1],'Vp'))
                    elif depthrange[d] <=50:
                        modvs2.append(crust.get_value(depthrange[d],lonP[-1],latP[-1],'Vs'))
                        modvp2.append(crust.get_value(depthrange[d],lonP[-1],latP[-1],'Vp'))
                    elif depthrange[d] <=2800:
                        dvs=mod.get_value(depthrange[d],lonPds[-1],latPds[-1],'dVs')
                        dvp=mod.get_value(depthrange[d],lonP[-1],latP[-1],'dVp')
                        vsref=taupmodel.model.s_mod.v_mod.evaluate_below(depthrange[d],'s')
                        vpref=taupmodel.model.s_mod.v_mod.evaluate_below(depthrange[d],'p')
                        modvs2.append((1+dvs/100.)*vsref)
                        modvp2.append((1+dvp/100.)*vpref)
 
                    else:
                        modvs2.append(taupmodel.model.s_mod.v_mod.evaluate_below(depthrange[d],'s'))
                        modvp2.append(taupmodel.model.s_mod.v_mod.evaluate_below(depthrange[d],'p'))
      
           
                    #Pwave
                    dtheta,dtime=godown(modvp2[d-1],modrad2[d-1],modvp2[d],modrad2[d],Pslow)
                    if np.isnan(dtheta):
                        dtheta=0.0
                        dtime=0.0
                    #print dtheta, dtime, modrad[d]
                    xp2.append(xp2[d-1]+dtheta)
                    point=line.Position(xp2[-1]*6371.e3)
                    #print slat,slon,point['lat2'],point['lon2']
                    latP.append(point['lat2'])
                    lonP.append(point['lon2'])
                    tpsum2.append(tpsum2[d-1]+dtime)

                    #Swave
                    dtheta,dtime=godown(modvs2[d-1],modrad2[d-1],modvs2[d],modrad2[d],Pdslow[d]) 
                    if np.isnan(dtheta):
                        dtheta=0.0
                        dtime=0.0
                    xs2.append(xs2[d-1]+dtheta)
                    point=line.Position(xs2[-1]*6371.e3)
                    latPds.append(point['lat2'])
                    lonPds.append(point['lon2'])
                    tssum2.append(tssum2[d-1]+dtime)

                    # calculate correction for distance between P and Pds where they hit depth d
                    r_pslow=Pslow*180./np.pi

                    dtcorr2.append((xp2[d]-xs2[d])*r_pslow)
                    tsum2.append(tssum2[d]-tpsum2[d]+dtcorr2[d])

                print(' time 3D is ', time.time()-start)
                #plt.plot(seis[0].conversions['ak135']['depths'],seis[0].conversions['ak135']['time'],'r-', linewidth=1)
                #plt.plot(depthrange,tsum2,'b-', linewidth=1)
                #plt.show()

                # save solutions to dictionary
                if not hasattr(seis[0],'conversions'):
                    seis[0].conversions=dict()

                seis[0].conversions['AFRP20CR']=dict()
                seis[0].conversions['AFRP20CR']['depths']=depthrange
                seis[0].conversions['AFRP20CR']['time']=np.array(tsum2)
                seis[0].conversions['AFRP20CR']['depthsfortime']=np.interp(seis[0].jgf1['time'],tsum2,depthrange)
                seis[0].conversions['AFRP20CR']['latitudes']=np.array(latPds)
                seis[0].conversions['AFRP20CR']['longitudes']=np.array(lonPds)
                seis[0].conversions['AFRP20CR']['latfortime']=np.interp(seis[0].jgf1['time'],tsum2,latPds)
                seis[0].conversions['AFRP20CR']['lonfortime']=np.interp(seis[0].jgf1['time'],tsum2,lonPds)
                    
                seis.write(stalist[i],format='PICKLE')
                print('AFRP20CR 3D depth conversion written for ' + stalist[i])

            else:
                print('AFRP20CR 3D depth conversion found for ' + stalist[i])

# Initialize 3D model for crust and mantle
init_model()

# Write a loop through data here

alldirs=sorted(glob.glob('../Data/*'))

for dir in alldirs:
    compute_onestation(dir, filter=rffilter)
    print('AFRP20CR Depth conversion finished for: '+str(dir))
    
