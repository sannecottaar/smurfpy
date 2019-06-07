# RF CREATION SCRIPT FOR RF ANALYSIS#################

# This script deconvolves RT seismogram components to produce RF based on user inputs
# requires module receiver_function.py in same folder program is run in. Adds computed RF to pre-existing PICKLE event file
#

# Import Modules
import obspy
from obspy import read
from obspy.core import Stream
from obspy.core import trace
import matplotlib.pyplot as plt
import os.path
import time
import glob
import shutil
import numpy as np
from obspy import UTCDateTime
import receiver_function as rf
import subprocess
from obspy.taup import TauPyModel


direc = '../Data'
flag = 'SV'  # rdial (SV) or tranverse (ST) RF to compute
filt = 'jgf1'  # define frequency range

if filt == 'jgf1':
    # uses gaussian pulses when making receiver functions
    filttype = 'gauss'
    # gaussian filter width (0.4 Hz - 2.5 sec)
    filterconst = 0.4
    # set filter 1 from 5 to 100 sec with 2.5 sec gaussians
    fmin = 0.01
    fmax = .2
if filt == 'jgf2':
    filttype = 'gauss'
    filterconst = 1.0
    # set filter 1 from 2 to 50 sec with 1 sec gaussians
    fmin = 0.02
    fmax = .5
if filt == 'jgf3':
    filttype = 'gauss'
    filterconst = .2
    # set filter 1 from 2 to 50 sec with 1 sec gaussians
    fmin = 0.02
    fmax = .1
if filt == 'tff1':
    filttype = 'gauss'
    filterconst = 1.4
    fmin = 0.01
    fmax = .7
if filt == 'tff2':
    filttype = 'gauss'
    filterconst = 1.2
    fmin = 0.01
    fmax = .6
if filt == 'tff3':
    filttype = 'gauss'
    filterconst = 1.0
    fmin = 0.01
    fmax = .5
if filt == 'tff4':
    filttype = 'gauss'
    filterconst = 0.8
    fmin = 0.01
    fmax = .4
if filt == 'tff5':
    filttype = 'gauss'
    filterconst = 0.6
    fmin = 0.01
    fmax = .3


# loop through events
stalist = glob.glob(direc + '/*/*.PICKLE')

c = 0
# Loop through data
for i in range(len(stalist)):  # range(cat.count()):
    print(i,stalist[i])

    # read in current event
    seis = read(stalist[i], format='PICKLE')

    # extract predicted P arrival time and epicentral dist
    Ptime = seis[0].stats.traveltimes['P']  # set P arrival time
    distdg = seis[0].stats['dist']

    # check to see if RF has already been computed for this event
    if hasattr(seis[0], filt):
        print('done with', stalist[i])
    else:

        Ptime = seis[0].stats.event.origins[0].time + Ptime
        vertical = seis.select(channel='*HZ')[0]
        Pref = vertical.slice(Ptime -25.,
            Ptime + 150.)  # Cut out P arrival on vertical

        Pref.time = Pref.times() - 25.

        if flag == 'SV':
            radial = seis.select(channel='*HR')[0]
            SVref = radial.slice(
                Ptime - 25.,
                Ptime + 150.)  # Cut out P arrival on radial
        if flag == 'SH':
            transverse = seis.select(channel='*HT')[0]
            SVref = transverse.slice(
                Ptime - 25.,
                Ptime + 150.)  # Cut out P arrival on transverse, still named SVref for simplicitiy

        # Resample data to 10 samples/sec (instead of 40 samples/sec). This
        # greatly speeds up all the computations
        SVref.resample(10.)
        SVref.time = SVref.times() - 25.
        # Filter and taper data
        Pref.filter(
            'bandpass',
            freqmin=fmin,
            freqmax=fmax,
            corners=2,
            zerophase=True)
        SVref.filter(
            'bandpass',
            freqmin=fmin,
            freqmax=fmax,
            corners=2,
            zerophase=True)
        seis[0].stats['minfreq'] = fmin
        seis[0].stats['maxfreq'] = fmax
        Pref.taper(max_percentage=0.05, type='cosine')
        SVref.taper(max_percentage=0.05, type='cosine')

        RF = trace.Trace()

        # Define dictionary for receiver function and filter details

        if flag == 'SV':

            if filt == 'jgf1':
                seis[0].jgf1 = dict()
                out = seis[0].jgf1

            elif filt == 'jgf2':
                seis[0].jgf2 = dict()
                out = seis[0].jgf2

            elif filt == 'jgf3':
                seis[0].jgf3 = dict()
                out = seis[0].jgf3

            elif filt == 'tff1':
                seis[0].tff1 = dict()
                out = seis[0].tff1

            elif filt == 'tff2':
                seis[0].tff2 = dict()
                out = seis[0].tff2

            elif filt == 'tff3':
                seis[0].tff3 = dict()
                out = seis[0].tff3

            elif filt == 'tff4':
                seis[0].tff4 = dict()
                out = seis[0].tff4

            elif filt == 'tff5':
                seis[0].tff5 = dict()
                out = seis[0].tff5

            else:
                print('Define a new dictionary to store results in')
        if flag == 'SH':
            if filt == 'rff1':
                seis[0].rfshf1 = dict()
                out = seis[0].rfshf1
            elif filt == 'rff2':
                seis[0].rfshf2 = dict()
                out = seis[0].rfshf2
            elif filt == 'rff3':
                seis[0].rfshf3 = dict()
                out = seis[0].rfshf3

        out['filter'] = filttype
        out['filterconst'] = filterconst
        out['timeshift'] = 25.

        # check that both components are not zero length otherwise script will
        # crash
        if len(SVref.data) == 0 or len(Pref.data) == 0:
            os.remove(stalist[i])
            continue

        # Run iterative deconvolution
        # this is were the actual work happens
        try:
            RF.data, fitid = rf.iterative_deconvolution(
            SVref.data, Pref.data, maxbumps=200, dt=Pref.stats['delta'], filt=filttype, fmax=filterconst, timeshift=25.)
        except:
            os.remove(stalist[i])            
            continue
        # Normalize and switch polarity if needed
        indm = np.argmax(np.abs(RF.data[200:300]))
        RF.data = RF.data / RF.data[200 + indm]

        out['iterativedeconvolution'] = RF.data
        out['iterativedeconvolution_fit'] = fitid
        out['time'] = SVref.time

        seis.write(stalist[i], format='PICKLE')
