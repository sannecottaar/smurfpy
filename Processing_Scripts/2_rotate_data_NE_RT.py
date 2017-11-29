#####DATA PROCESSING FOR RF ANALYSIS#################
# This script: merges cut up components into one stream
#              trims components so they are all one length
#              downsamples data to 10sps
#              rotates seismograms from ZNE to ZRT orientations
#              Renames based on BAZ and epicentral distance infomation

# Outputs: -python stream objects in PICKLE format with a dictionary of header info for each event in a new folder leaving a copy of the un processed original data for future use


import obspy
from obspy import read
from obspy.core import Stream
from obspy.core import event
from obspy.taup.taup import getTravelTimes
from obspy import UTCDateTime
import obspy.signal
import obspy.signal.rotate
import matplotlib.pyplot as plt
import os.path
import time
import glob
import shutil
import numpy as np
import scipy
import sys,os

# Find list of stations directories
stations = glob.glob('DataRF/*')

# Loop through stations
for stadir in stations:
        print("STATION:", stadir)
        stalist=glob.glob(stadir+'/Originals/*PICKLE') 
        # make directory for processed data
        direc= stadir+'/Processed_originals'
        if not os.path.exists(direc):
                os.makedirs(direc)
        # Loop through events
        for s in range(len(stalist)):
            print(stalist[s])
            
            #try:
            onestation = read(stalist[s],format='PICKLE')

            #note stats of stream as merge gets rid of them
            BAZ=onestation[0].stats['baz']
            EVLA=onestation[0].stats['evla']
            EVLO= onestation[0].stats['evlo']
            EVDP= onestation[0].stats['evdp']
            STLA= onestation[0].stats['stla']
            STLO= onestation[0].stats['stlo']
            DIST= onestation[0].stats['dist']
            AZ= onestation[0].stats['az']
            STATION= onestation[0].stats['station']
            NETWORK= onestation[0].stats['network']
            EVENT= onestation[0].stats['event']
                
            #detrending
            onestation.detrend('linear')
        
            #merge gappy waveforms and overlap
            onestation.merge()

            #find streams with min/max start/end times
            start_cut=UTCDateTime(onestation[0].stats.starttime)
            end_cut=UTCDateTime(onestation[0].stats.endtime)
            for s1 in onestation:
                    print(s1.stats.sampling_rate)

                    start_time_stream=UTCDateTime(s1.stats.starttime)
                    end_time_stream=UTCDateTime(s1.stats.endtime)

                    print("STREAM:", s1)
                    print(start_time_stream)
                    print(end_time_stream)
                    if  start_time_stream > start_cut:

                         start_cut=start_time_stream
                    if  end_time_stream < end_cut:
                         end_cut=end_time_stream


            # Cut components to the same length
            onestation.trim(starttime=start_cut, endtime=end_cut)

            # Find if component names are differenst
            seisZ=onestation.select(channel='BHZ')           
            if len(seisZ)>1:
                seisZ=seisZ.select(sampling_rate=20.0)
            if len(seisZ)>1:
                seisZ=seisZ.select(location='10')                

            seisN=onestation.select(channel='BHN')
            
            if len(seisN)==0:
                seisN=onestation.select(channel='BH2')
            if len(seisN)>1:
                seisN=seisN.select(sampling_rate=20.0)                 
            if len(seisN)>1:
                seisN=seisN.select(location='10')                
            seisE=onestation.select(channel='BHE')        
            if len(seisE)==0:
                seisE=onestation.select(channel='BH1')
            if len(seisE)>1:
                seisE=seisE.select(sampling_rate=20.0)                
            if len(seisE)>1:
                seisE=seisE.select(location='10')


            # rotate components to from North and East to Radial and Transverse
            [seisRtmp,seisTtmp]=obspy.signal.rotate.rotate_ne_rt(seisN[0].data,seisE[0].data,BAZ)
            seisR=seisN[0].copy()
            seisR.stats['channel']='BHR'
            seisR.data=seisRtmp
            seisT=seisN[0].copy()
            seisT.stats['channel']='BHT'
            seisT.data=seisTtmp

  
            
            # produce new stream with Vertical Radial and Transverse
            seisnew=Stream()
            seisnew.append(seisZ[0])
            seisnew.append(seisR)
            seisnew.append(seisT)
            for s in seisnew:
                print(s.stats.npts,s.stats.delta)   
                #resample to 10 samples/s
                s.resample(10)
                print("resampled to:",s.stats.npts,s.stats.delta)
                
            # Copy values into stats 
            seisnew[0].stats['evla'] =EVLA
            seisnew[0].stats['evlo'] =EVLO
            seisnew[0].stats['evdp'] =EVDP
            seisnew[0].stats['stla'] =STLA
            seisnew[0].stats['stlo'] =STLO
            seisnew[0].stats['dist'] =DIST
            seisnew[0].stats['az']   =AZ
            seisnew[0].stats['baz']  = BAZ
            seisnew[0].stats['station'] =STATION
            seisnew[0].stats['network'] =NETWORK
            seisnew[0].stats['event'] =EVENT

            ## write out in pickle format (with BAZ and DIST in name)
            baz_str=format(BAZ,'.3f').zfill(7)
            dist_str=format(DIST,'.3f').zfill(7)
            ev_str=str(seisZ[0].stats.starttime)
            filename=stadir+'/'+baz_str+'_'+dist_str+'_'+str(seisZ[0].stats.starttime)+'.PICKLE'

            seisnew.write(filename,'PICKLE')
            shutil.move(stalist[i], direc)
            #except:
            #    print('FAILED for', stalist[s])






