# DATA DOWNLOAD FOR RF ANALYSIS#################
# This script selects appropriate events and stations within a set area based on user inputs
# Data is downloaded from the IRIS server

# Inputs: -search area lat/lon
#         -start and end times (from command line if commented in)
#         -epicentral dist of ev/station
#         -event magnitude range
#         -trace lenth
#         -data filter band
#         -station networks to search

# Outputs: -python stream objects in PICKLE format with a dictionary of
# header info for each event

# NOTE: Recommended to run in yearly blocks max 3 in parallel, otherwise
# the event catalogues get so big that the connection with IRIS times out.
# Also note that it is important to check IRIS stations and pick out
# relevant networks otherwise program wastes a lot of time look in single
# component or infra sound only networks. MUST make folder 'DataRF' in
# current working directory program will be run in - the structure will be
# used in all following scripts


# Usage is: python download_data_per_station.py

import obspy
from obspy.clients.fdsn import Client as IRISClient
from obspy import UTCDateTime
import os.path
import obspy.geodetics.base
import numpy as np
import obspy.geodetics
import sys

# Command line help
if len(sys.argv) > 1:
    print('\n')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print(sys.argv[0])
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('Description:           select and download appropriate events and stations based on user inputs')
    print('Inputs:                search area lat/lon, start and end times (in function), epicentral dist of ev/station,')
    print('                       event magnitude range, trace length, data filter band, station networks to search, dataclient')
    print('Outputs:               python stream objects in PICKLE format with a dictionary of header info for each event.')
    print('                       Saves to ../Data/NT.STA/Originals/\n')
    print('Usage:                 >> python3 1_download_data_per_station.py')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('\n')
    sys.exit()

def download_data(start, end):
    # Load IRIS client. If using another client (e.g. RESIF, GEOFON...) edit dataclient 
    # accordingly
    irisclient = IRISClient("IRIS")
    dataclient = IRISClient("IRIS")
    
    # Station parameters
    # longitude and latitude bounds for stations (currently set  to Samoa)
    lonmin = - 175  
    lonmax = -165
    latmin = -15
    latmax = -10

    # Set time bands (at the bottom of the script)
    starttime = UTCDateTime(start)
    endtime = UTCDateTime(end)

    # Event parameters
    radmin = 30  # Minimum radius for teleseismic earthquakes
    radmax = 90  # Maximum radius (further distances interact with the core-mantle boundary
    minmag = 5.5  # Minumum magnitude of quake
    maxmag = 8.0  # Maximum magnitude of quake
    lengthoftrace = 30. * 60.  # 30 min

    # Define a filter band to prevent amplifying noise during the deconvolution
    # not being used, as responses are not being removed
    fl1 = 0.005
    fl2 = 0.01
    fl3 = 2.
    fl4 = 20.
    # Find all suitable stations on IRIS
    inventory = irisclient.get_stations(
        minlatitude=latmin,
        maxlatitude=latmax,
     minlongitude=lonmin,
     maxlongitude=lonmax,
     starttime=starttime,
     endtime=endtime,
     channel='*H*',
     level="channel")
    count = 0
    stacount = 0
    for nw in inventory:
        for sta in nw:
            name = nw.code + '.' + sta.code
            print(sta)

            # Make directories for data
            direc = '../Data/' + name
            if not os.path.exists(direc):
                os.makedirs(direc)
            direc = direc + '/Originals'
            if not os.path.exists(direc):
                os.makedirs(direc)

            # Find events
            print(sta.code, nw.code)
            mintime = sta.start_date
            if mintime < starttime:
                mintime = starttime
            maxtime = sta.end_date
            if maxtime > endtime:
                maxtime = endtime
            # Find all suitable events
            cat = irisclient.get_events(
                latitude=sta.latitude,
                longitude=sta.longitude,
                minradius=radmin,
                maxradius=radmax,
                starttime=mintime,
                endtime=maxtime,
                minmagnitude=minmag)
            print(len(cat))
            for ev in cat:
                evtime = ev.origins[0].time
                t = UTCDateTime(evtime)
                seis = []
                try:
                    print(
                        'downloading',
                        nw.code,
                        sta.code,
                        "*",
                        "BH*",
                        evtime,
                        evtime +
                        lengthoftrace)
                    seis = dataclient.get_waveforms(
                        nw.code,
                        sta.code,
                        "*",
                        "BH*",
                        evtime,
                        evtime + lengthoftrace,
                        attach_response=True)
                except:
                    try:
                        print(
                            'Attempting',
                            nw.code,
                            sta.code,
                            "*",
                            "HH*",
                            evtime,
                            evtime +
                            lengthoftrace)
                        seis = dataclient.get_waveforms(
                            nw.code,
                            sta.code,
                            "*",
                            "HH*",
                            evtime,
                            evtime + lengthoftrace,
                            attach_response=True)
                    except:
                        print('failed',nw.code,sta.code, "*","*H*",evtime,evtime +
                            lengthoftrace)
                
                if len(seis) > 1:
                    evtlatitude = ev.origins[0]['latitude']
                    evtlongitude = ev.origins[0]['longitude']
                    try:
                        evtdepth = ev.origins[0][
                            'depth'] / 1.e3  # convert to km from m
                    except:
                        print('failed to get true depth')
                        evtdepth = 30.  # This is a bit of a hack, might need better solution

                    # Compute distances azimuth and backazimuth
                    distm, az, baz = obspy.geodetics.base.gps2dist_azimuth(
                        evtlatitude, evtlongitude, sta.latitude, sta.longitude)
                    distdg = distm / (6371.e3 * np.pi / 180.)

                    # Remove station instrument response. + Output seismograms
                    # in displacement and filter with corner frequences fl2,
                    # fl3
                    try:
                        seis.remove_response(
                            output='DISP', pre_filt=[fl1, fl2, fl3, fl4])
                    except:
                        print('failed to remove response')

                    # Put in various event and station characteristics into the
                    # 'stats'-dictionairy
                    seis[0].stats['evla'] = evtlatitude
                    seis[0].stats['evlo'] = evtlongitude
                    seis[0].stats['evdp'] = evtdepth
                    seis[0].stats['stla'] = sta.latitude
                    seis[0].stats['stlo'] = sta.longitude
                    seis[0].stats['dist'] = distdg
                    seis[0].stats['az'] = az
                    seis[0].stats['baz'] = baz
                    seis[0].stats['station'] = sta.code
                    seis[0].stats['network'] = nw.code
                    seis[0].stats['event'] = ev
                    
                    # Append azimuth and dip to each channels stats directory   
                    for cha in range(len(seis)):
                        location = seis[cha].stats['location']
                        channel = seis[cha].stats['channel']
                        identifier = str(nw.code) + "." + str(sta.code) + "." \
                        + str(location) + "." + str(channel)                       
                        orientation = nw.get_orientation(str(identifier), t)
                        seis[cha].stats['orientation'] = orientation['azimuth']
                        seis[cha].stats['dip'] = orientation['dip']
                        # Add the station elevation to the vertical component
                        if cha == 0:
                            inv = inventory.select(network=str(nw.code), station=str(sta.code))
                            stel=inv.get_coordinates(identifier,t)['elevation']/1000 # stel in kilometers, positive is upward topography.
                            seis[cha].stats['stel'] = stel


                    # Write out to file
                    filename = direc + '/' + \
                        str(seis[0].stats.starttime) + '.PICKLE'
                    print('writing to ', filename)
                    count = count + 1
                    print(count)
                    seis.write(filename, format='PICKLE')

    print(starttime,endtime,str(count),'seismograms found for',str(stacount),'stations')


download_data('1992-01-01', '2018-01-01')
