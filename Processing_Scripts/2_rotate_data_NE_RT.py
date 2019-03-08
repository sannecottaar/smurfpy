# DATA PROCESSING FOR RF ANALYSIS#################
# This script: merges cut up components into one stream
#              trims components so they are all one length
#              downsamples data to 10sps
#              rotates seismograms from ZNE to ZRT orientations
#              Renames based on BAZ and epicentral distance infomation

# Outputs: -python stream objects in PICKLE format with a dictionary of
# header info for each event in a new folder leaving a copy of the un
# processed original data for future use


import obspy
from obspy import read
from obspy.core import Stream
from obspy import UTCDateTime
import obspy.signal
import obspy.signal.rotate
import os.path
import glob
import shutil
import os

no_trace = 0
no_z = 0

# Find list of stations directories
stations = glob.glob('../Data/*')

f = open("error_file.txt", 'w')
all_data = 0
# Loop through stations
for stadir in stations:
    print("STATION:", stadir)
    stalist = glob.glob(stadir + '/Originals/*PICKLE')
    # make directory for processed data
    direc = stadir + '/Processed_originals'
    direc1 = stadir + '/No_Z_component'
    direc2 = stadir + '/No_NE_component'
    if not os.path.exists(direc):
        os.makedirs(direc)
    if not os.path.exists(direc1):
        os.makedirs(direc1)
    if not os.path.exists(direc2):
        os.makedirs(direc2)
    # Loop through events
    for s in range(len(stalist)):
        print(s,stalist[s])
        
        all_data += 1

        onestation = read(stalist[s], format='PICKLE')

        # note stats of stream as merge gets rid of them
        BAZ = onestation[0].stats['baz']
        EVLA = onestation[0].stats['evla']
        EVLO = onestation[0].stats['evlo']
        EVDP = onestation[0].stats['evdp']
        STLA = onestation[0].stats['stla']
        STLO = onestation[0].stats['stlo']
        DIST = onestation[0].stats['dist']
        AZ = onestation[0].stats['az']
        STATION = onestation[0].stats['station']
        NETWORK = onestation[0].stats['network']
        EVENT = onestation[0].stats['event']

        # detrending
        onestation.detrend('linear')

        # merge gappy waveforms and overlap
        onestation.merge()

        # find streams with min/max start/end times
        start_cut = UTCDateTime(onestation[0].stats.starttime)
        end_cut = UTCDateTime(onestation[0].stats.endtime)
        for s1 in onestation:
            print(s1.stats.sampling_rate)

            start_time_stream = UTCDateTime(s1.stats.starttime)
            end_time_stream = UTCDateTime(s1.stats.endtime)

            print("STREAM:", s1)
            #print(start_time_stream)
            #print(end_time_stream)

            #print(start_cut)
            #print(end_cut)

            if start_time_stream > start_cut:

                start_cut = start_time_stream
            if end_time_stream < end_cut:
                end_cut = end_time_stream

            #print(start_cut)
            #print(end_cut)

        # Cut components to the same length
        onestation.trim(starttime=start_cut, endtime=end_cut)
        
        # Find component names
        seisZ = onestation.select(channel='BHZ')
        seisN = onestation.select(channel='BHN')
        seisE = onestation.select(channel='BHE')
        
        N_E = True
        
        # Get data if BH1/2 rather than BHN/E
        if len(seisN) == 0 or len(seisE) == 0:
            N_E = False
            seisN = onestation.select(channel='BH2')
            seisE = onestation.select(channel='BH1')
        
        # Check North and East exist
        if len(seisN) >=1 and len(seisE)>=1:
            
            # Check Z exists
            if len(seisZ) >=1:
                
                # Make sure there are no more than 3 channels
                if len(seisZ) > 1 or len(seisN) > 1 or len(seisE) > 1:
                    if seisZ[0].stats['location'] == '00':
                        seisZ = seisZ.select(location='00')
                        seisN = seisN.select(location='00')
                        seisE = seisE.select(location='00')
                        other_loc1 = '10'
                        other_loc2 = '01'
                    elif seisZ[0].stats['location'] == '10':
                        seisZ = seisZ.select(location='10')
                        seisN = seisN.select(location='10')
                        seisE = seisE.select(location='10')
                        other_loc1 = '01'
                        other_loc2 = '00'
                    elif seisZ[0].stats['location'] == '01':
                        seisZ = seisZ.select(location='01')
                        seisN = seisN.select(location='01')
                        seisE = seisE.select(location='01')
                        other_loc1 = '00'
                        other_loc2 = '10'
                
                # Check thoroughly to get 3 components
                # Only check if N & E are both > 1
                if N_E:
                    if len(seisZ) == 0 or len(seisN) == 0 or len(seisE) == 0:
                        seisZ = onestation.select(channel='BHZ',location=other_loc1)
                        seisN = onestation.select(channel='BHN',location=other_loc1)
                        seisE = onestation.select(channel='BHE',location=other_loc1)
                        
                    if len(seisZ) == 0 or len(seisN) == 0 or len(seisE) == 0:
                        seisZ = onestation.select(channel='BHZ',location=other_loc2)
                        seisN = onestation.select(channel='BHN',location=other_loc2)
                        seisE = onestation.select(channel='BHE',location=other_loc2)
                    
                # Always check this in case BHE & BHN have different location    
                if len(seisZ) == 0 or len(seisN) == 0 or len(seisE) == 0:
                    seisZ = onestation.select(channel='BHZ',location=other_loc1)
                    seisN = onestation.select(channel='BH2',location=other_loc1)
                    seisE = onestation.select(channel='BH1',location=other_loc1)
                    
                if len(seisZ) == 0 or len(seisN) == 0 or len(seisE) == 0:
                    seisZ = onestation.select(channel='BHZ',location=other_loc2)
                    seisN = onestation.select(channel='BH2',location=other_loc2)
                    seisE = onestation.select(channel='BH1',location=other_loc2)
                    
                
                # Make sure components are of the same length
                while len(seisN[0].data) != len(seisE[0].data):
                    if len(seisN[0].data) > len(seisE[0].data):
                        seisN[0].data = seisN[0].data[:-1]
                    if len(seisN[0].data) < len(seisE[0].data):
                        seisE[0].data = seisE[0].data[:-1]
                    
                # Make sure North component is North and East component is East
                if 45 < seisN[0].stats['orientation'] < 135 or 225 < seisN[0].stats['orientation'] < 315:
                    seisN1 = seisE
                    seisE1 = seisN
                    seisN = seisN1
                    seisE = seisE1
            
                
            
                #Rotate BH1/BH2 to BHN/BHE
                if not N_E:
                    
                    orientationN = seisN[0].stats['orientation']
                    dipN = seisN[0].stats['dip']
            
                    orientationE = seisE[0].stats['orientation']
                    dipE = seisE[0].stats['dip']
                
                    orientationZ = seisZ[0].stats['orientation']
                    dipZ = seisZ[0].stats['dip']
                    
                    print('Rotating to N/E')
                    [seisZtmp, seisNtmp, seisEtmp] = obspy.signal.rotate.rotate2zne(
                        seisZ[0].data, orientationZ, dipZ,
                        seisN[0].data, orientationN, dipN,
                        seisE[0].data, orientationE, dipE)
                    seisN = seisN[0].copy() 
                    seisN.stats['channel'] = 'BHN'
                    seisN.data = seisNtmp
                    seisE = seisE[0].copy()
                    seisE.stats['channel'] = 'BHE'
                    seisE.data = seisEtmp
                    
                # rotate components to from North and East to Radial and Transverse
        
                    [seisRtmp, seisTtmp] = obspy.signal.rotate.rotate_ne_rt(
                        seisN.data, seisE.data, BAZ)
                    seisR = seisN.copy()
                    seisT = seisN.copy()
                  
                else:
                    [seisRtmp, seisTtmp] = obspy.signal.rotate.rotate_ne_rt(
                        seisN[0].data, seisE[0].data, BAZ)
                    seisR = seisN[0].copy()
                    seisT = seisN[0].copy()
                    
                seisR.stats['channel'] = '*HR'
                seisR.data = seisRtmp
                seisT.stats['channel'] = '*HT'
                seisT.data = seisTtmp
        
                # produce new stream with Vertical Radial and Transverse
                seisnew = Stream()
                seisnew.append(seisZ[0])
                seisnew.append(seisR)
                seisnew.append(seisT)
                
                for st in seisnew:
                    print(st.stats.npts, st.stats.delta)
                # resample to 10 samples/s
                    st.resample(10)
                    print("resampled to:", st.stats.npts, st.stats.delta)
        
            # Copy values into stats
                seisnew[0].stats['evla'] = EVLA
                seisnew[0].stats['evlo'] = EVLO
                seisnew[0].stats['evdp'] = EVDP
                seisnew[0].stats['stla'] = STLA
                seisnew[0].stats['stlo'] = STLO
                seisnew[0].stats['dist'] = DIST
                seisnew[0].stats['az'] = AZ
                seisnew[0].stats['baz'] = BAZ
                seisnew[0].stats['station'] = STATION
                seisnew[0].stats['network'] = NETWORK
                seisnew[0].stats['event'] = EVENT
        
            # write out in pickle format (with BAZ and DIST in name)
                baz_str = '%03d' %(BAZ)
                dist_str = '%02d' %(DIST)
        
                ev_str = str(seisZ[0].stats.starttime)
                filename = stadir + '/' + baz_str + '_' + dist_str + \
                    '_' + str(seisZ[0].stats.starttime) + '.PICKLE'
        
                seisnew.write(filename, 'PICKLE')
                shutil.move(stalist[s], direc)
                
            else:
                #print('No Z')
                no_z += 1
                f.write("Station: " + str(stadir) + ", Event: " + str(s) + '\n')
                f.write("Failed on Z component" + '\n')
                shutil.move(stalist[s], direc1)
    

        else:
            #print('Trace error in seisN/seisE')
            no_trace += 1
            f.write("Station: " + str(stadir) + ", Event: " + str(s) + '\n')
            f.write("Failed on N/E component" + '\n')
            shutil.move(stalist[s], direc2)

print(all_data)
print(no_trace)
print(no_z)

f.close()

