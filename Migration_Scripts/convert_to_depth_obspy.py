############## Convert RF from time to depth using 1D model ##############
# Converts RF from time to depth using obspy and taup
# 1D model defined in 'taupmodel' below imports
# Change station directory at base of script to point to your station dir
# Outputs dictionary seis[0].conversions['<nameof1Dmodel>']
##########################################################################

import sys,os 
sys.path.append('/raid2/sc845/Python/geographiclib-1.34/')
from obspy import read
from obspy.core import trace
import os.path
import glob
import numpy as np
import subprocess
from obspy.taup import TauPyModel


taupmodel = TauPyModel(model='./prem_added_discon_taup.npz')

def compute_onestation(dr,filter='jgf1'):

    # Loop through receiver functions in station directory
    stalist=glob.glob(dr+'/*.PICKLE')

    # Loop through list
    for i in range(len(stalist)): 
        print(i, stalist[i])
        seis=read(stalist[i],format='PICKLE') # Read in receiver function
        RFtrace = getattr(seis[0],filter)
        if not 'depths' in RFtrace:
          if seis[0].stats['dist']> 30.  and seis[0].stats['dist'] <90.:
            print('converting for ', stalist[i])
            dist=seis[0].stats['dist']                     # Set distance
            depth= seis[0].stats['evdp'] # set depth
            slat=float(seis[0].stats['stla'])
            slon=float(seis[0].stats['stlo'])
            elat= seis[0].stats['event'].origins[0]['latitude']
            elon= seis[0].stats['event'].origins[0]['longitude']

            RF=trace.Trace()
            RF.data= RFtrace['iterativedeconvolution']      # Define iterative deconvolution as a trace

   
            # Initialize
            refdepths=[]
            reftimes=[]
            reftakeoff=[]
            refincident=[]
            refslow=[]

            # Get slownesses and incident angles based on PREM using Taup in Obspy
            phases =  ['P','Pms','P220s','P310s','P400s','P550s','P670s','P971s','P1171s']
            phdepths = [0,24.4,220,310,400,550,670,971,1171]
            for p,ph in enumerate(phases):
                arr = taupmodel.get_travel_times(depth,dist,[ph])
                if p ==0:
                    Preftime=arr[0].time
                    Prefslow=arr[0].ray_param
                    Preftakeoff= arr[0].takeoff_angle
                    Prefincident = arr[0].incident_angle
                refdepths.append(phdepths[p])
                reftimes.append(arr[0].time-Preftime)
                reftakeoff.append(arr[0].takeoff_angle-Preftakeoff)
                refincident.append(arr[0].incident_angle-Prefincident)
                refslow.append((arr[0].ray_param-Prefslow)*np.pi/180.)
            depthspace=np.linspace(-100,1200,1300)
            # Fit polynomials
            takeoff_poly=np.polyfit(refdepths[1:],reftakeoff[1:],3)
            takeoff=np.poly1d(takeoff_poly)
            slow_poly=np.polyfit(refdepths,refslow,4)
            slowness=np.poly1d(slow_poly)
            time_poly=np.polyfit(refdepths,reftimes,4)
            timepl=np.poly1d(time_poly)

            # Save solutions to dictionary
            if not hasattr(seis[0],'conversions'):
                seis[0].conversions=dict()
            seis[0].conversions['prem']=dict()
            seis[0].conversions['prem']['depths']=depthspace
            seis[0].conversions['prem']['takeoff']=takeoff(depthspace)
            seis[0].conversions['prem']['slowness']=slowness(depthspace)
            seis[0].conversions['prem']['time']=timepl(depthspace)




            # Interpolate latitude and longitudes  at the discontinuities
            # Unfortunately obpsy taup cannot do this...yet. So this will be slow
            discon=['P','P220s','P400s','P670s', 'P1171s']
            x=[0,220.,400.,670., 1171.]

            xtime=[0]
            xlat=[slat]
            xlon=[slon]
            for d in range(1,len(discon)):
  
                taup_660=['taup_pierce -mod prem_added_discon_taup.nd -h '+ str(depth) + ' -sta '+ str(slat)+' ' +str(slon)+ ' -evt '+ str(elat) + ' ' + str(elon) +' -ph '+discon[d] + ' -Pierce ' + str(x[d])+',0 -nodiscon']
                out=subprocess.check_output(taup_660,shell=True)
                t=out.split()
                xtime.append(float(t[-3])-float(t[-8]))
                xlat.append(float(t[-7]))
                xlon.append(float(t[-6]))

            traveltime_poly=np.polyfit(x,xtime,4)
            traveltime=np.poly1d(traveltime_poly)
            lats_poly=np.polyfit(x,xlat,4)
            lats=np.poly1d(lats_poly)
            lons_poly=np.polyfit(x,xlon,4)
            lons=np.poly1d(lons_poly)


            # Convert this trace to depths and save
            depths=np.interp(RFtrace['time'],seis[0].conversions['prem']['time'],seis[0].conversions['prem']['depths'])
            latitudes=lats(depthspace)
            longitudes=lons(depthspace)

            seis[0].conversions['prem']['depthsfortime']=depths
            seis[0].conversions['prem']['latitudes']=latitudes
            seis[0].conversions['prem']['longitudes']=longitudes
            seis[0].conversions['prem']['latfortime']=np.interp(RFtrace['time'],seis[0].conversions['prem']['time'], latitudes)
            seis[0].conversions['prem']['lonfortime']=np.interp(RFtrace['time'],seis[0].conversions['prem']['time'],longitudes)

            seis.write(stalist[i],format='PICKLE')
          else:
            os.remove(stalist[i])               

stations=glob.glob('../Data/*')

for dr in stations:
    compute_onestation(dr)
