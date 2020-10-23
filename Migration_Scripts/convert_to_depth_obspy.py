############## Convert RF from time to depth using 1D model ##############
# Converts RF from time to depth using obspy and taup
# 1D model defined in 'taupmodel' below imports
# Change station directory at base of script to point to your station dir
# Outputs dictionary seis[0].conversions['<nameof1Dmodel>']
##########################################################################

import sys,os 
from obspy import read
from obspy.core import trace
import os.path
import glob
import numpy as np
import subprocess
from obspy.taup import TauPyModel

# Command line help
if len(sys.argv) != 2 or str(sys.argv[1]).lower() == 'help':
    print('\n')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print(sys.argv[0])
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('Description:           Convert RF from time to depth using 1D model (coded for Prem)')
    print('Inputs:                Filter band, 1D velocity model')
    print("Outputs:               Adds dictionary seis[0].conversions['<nameof1Dmodel>'] to each Pickle file\n")
    print('Usage:                 >> python3 convert_to_depth_obspy.py filterband')
    print('Options [1]:           jgf1, jgf2, jgf3, tff1, tff2, tff3, tff4 or tff5 [str]')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('\n')
    sys.exit()

PREM=True
ak135=False
rffilter=str(sys.argv[1])

if PREM:
    mod_1D='prem'
    taupmodel = TauPyModel(model='../Tools/MODELS/PREM_FILES/prem_added_discon_taup.npz')
elif ak135:
    mod_1D='ak135'
    #obspy.taup.taup_create.build_taup_model('./ak135_added_discon_taup.tvel','./ak135_added_discon_taup.npz')
    taupmodel = TauPyModel(model='../Tools/MODELS/ak135_FILES/ak135_added_discon_taup.npz')


def compute_onestation(dr,filter='jgf1'):

    # Loop through receiver functions in station directory
    station_file = open(dr + '/selected_RFs_'+str(filter)+'.dat','r')
    stalist=station_file.read().strip().split()

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

            
            if PREM:
                # Get slownesses and incident angles based on PREM using Taup in Obspy
                phases =  ['P','Pms','P220s','P310s','P400s','P550s','P670s','P971s','P1171s']
                phdepths = [0,24.4,220,310,400,550,670,971,1171]
            if ak135:
                # get slownesses and incident angles based on ak135 using Taup in Obspy
                phases =  ['P','Pms','P210s','P410s','P660s','P1156s']
                phdepths = [0,35,210,410,660,1156]


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
            seis[0].conversions[str(mod_1D)]=dict()
            seis[0].conversions[str(mod_1D)]['depths']=depthspace
            seis[0].conversions[str(mod_1D)]['takeoff']=takeoff(depthspace)
            seis[0].conversions[str(mod_1D)]['slowness']=slowness(depthspace)
            seis[0].conversions[str(mod_1D)]['time']=timepl(depthspace)



            if PREM:
                # Interpolate latitude and longitudes  at the discontinuities
                # Unfortunately obpsy taup cannot do this...yet. So this will be slow
                discon=['P','P220s','P400s','P670s', 'P1171s']
                x=[0,220.,400.,670., 1171.]
            if ak135:
                discon=['P','P210s','P410s','P660s','P1156s']
                x=[0,210.,410.,660.,1156.]

            xtime=[0]
            xlat=[slat]
            xlon=[slon]
            for d in range(1,len(discon)):
                if PREM:
                    taup_660=['taup_pierce -mod ../Tools/MODELS/PREM_FILES/prem_added_discon_taup.nd -h '+ str(depth) + ' -sta '+ str(slat)+' ' +str(slon)+ ' -evt '+ str(elat) + ' ' + str(elon) +' -ph '+discon[d] + ' -Pierce ' + str(x[d])+',0 -nodiscon']
                if ak135:
                    taup_660=['taup_pierce -mod ../Tools/MODELS/ak135_FILES/ak135_added_discon_taup.tvel -h '+ str(depth) + ' -sta '+ str(slat)+' ' +str(slon)+ ' -evt '+ str(elat) + ' ' + str(elon) +' -ph '+discon[d] + ' -Pierce ' + str(x[d])+',0 -nodiscon']
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
            depths=np.interp(RFtrace['time'],seis[0].conversions[str(mod_1D)]['time'],seis[0].conversions[str(mod_1D)]['depths'])
            latitudes=lats(depthspace)
            longitudes=lons(depthspace)

            seis[0].conversions[str(mod_1D)]['depthsfortime']=depths
            seis[0].conversions[str(mod_1D)]['latitudes']=latitudes
            seis[0].conversions[str(mod_1D)]['longitudes']=longitudes
            seis[0].conversions[str(mod_1D)]['latfortime']=np.interp(RFtrace['time'],seis[0].conversions[str(mod_1D)]['time'], latitudes)
            seis[0].conversions[str(mod_1D)]['lonfortime']=np.interp(RFtrace['time'],seis[0].conversions[str(mod_1D)]['time'],longitudes)

            seis.write(stalist[i],format='PICKLE')
          else:
            os.remove(stalist[i])               

stations=glob.glob('../Data/*')

for dr in stations:
    compute_onestation(dr,filter=rffilter)
