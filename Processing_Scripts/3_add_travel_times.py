#####TRAVEL TIME PREDICTION SCRIPT FOR RF ANALYSIS#################
# This script computes predicted travel times for user defined phases based on TauP
# predicted times are added to the events header information

# Example: python3 3_add_travel_times P S P660s P410s#####DATA PROCESSING FOR RF ANALYSIS#################


import obspy
from obspy import read
from obspy.core import Stream
from obspy.core import event
from obspy import UTCDateTime 
import datetime
import time
import obspy.signal
import obspy.signal.rotate
import os.path
import glob
import shutil
import numpy as np
import scipy #what does this do?
import sys,os
import obspy
from obspy import read
from obspy.core import Stream
from obspy.core import event
from obspy.taup.taup import getTravelTimes
from obspy import UTCDateTime
import obspy.signal
import matplotlib.pyplot as plt
import os.path
import time
import glob
import shutil
import numpy as np
import scipy
from obspy.io.xseed import Parser
from obspy.clients.arclink import Client as ARCLINKClient
from obspy.clients.fdsn import Client as IRISClient
from subprocess import call
import subprocess
from obspy.taup.taup import getTravelTimes
import sys
from obspy.taup import TauPyModel
model=TauPyModel(model="prem")

# Find list of stations directories
stations = glob.glob('DataRF/*')

# Add phase names as additional arguments (these are the TauP phase names)
phase=[]
for i in range(1,len(sys.argv)):
    phase.append(sys.argv[i])
count=0
# Loop through stations
for stadir in stations:
        print("STATION:", stadir)
        stalist=glob.glob(stadir+'/*PICKLE') 
        # make directory for processed data
        direc= stadir+'/Travel_time_added'

        if not os.path.exists(direc):
                os.makedirs(direc)
        # Loop through events
        for sta in stalist:
            print(sta)
            seis = read(sta,format='PICKLE')
            
            try:
                    if not hasattr(seis[0].stats,'traveltimes'):
                            seis[0].stats.traveltimes=dict()
                            #print('making dict')
                            #Loop through phases and call TauP_time to get traveltime
                    for ph in range(len(phase)):
                        print(phase[ph])
                        if not phase[ph] in seis[0].stats.traveltimes:
                                # extract only first time value and if value does not exists ttime=None
                                try:
                                        evdep_km=seis[0].stats['evdp']
                                        #print(evdep_km)
                                        dist_deg=seis[0].stats['dist']
                                        #print(dist_deg)
                                        arrivals = model.get_travel_times(source_depth_in_km=evdep_km,distance_in_degree=dist_deg,phase_list=[phase[ph]])
                                        #print(arrivals)
                                        ttime=arrivals[0].time
                                        print(phase[ph],ttime)
                                except:
                                        ttime=None
                                        print(phase[ph],ttime)
                                seis[0].stats.traveltimes[phase[ph]]=ttime
                        else:
                                print('already calculated ttime for', phase[ph],'\t', s)# Write out seismogram again ## write out in pickle format (with BAZ and DIST in name)
                    BAZ=seis[0].stats['baz']
                    DIST=seis[0].stats['dist']
                    print(BAZ,DIST)
                    baz_str=format(BAZ,'.3f').zfill(7)
                    dist_str=format(DIST,'.3f').zfill(7)
                    ev_str=str(seis[0].stats.starttime)
                    #filename=stadir+'/'+baz_str+'_'+dist_str+'_'+ev_str+'.PICKLE'
                    seis.write(sta,'PICKLE')
                    count += 1
                    print('Wrote to File', sta, ' number', count)
                    shutil.copy(sta, direc)
                    print('Copied to:', direc)
                    os.remove(sta)
                    print('Now deleting ', sta)
            except:
                    print('FAILED for', sta)
                    print(sys.exc_info)





