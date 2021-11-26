import matplotlib.pyplot as plt
import numpy as np
import glob
import os

# loop thru stations
#   for a given station save its name and find length of rf .dat file
# plot lengths

def rfsperstation(Data, noisefilter, filt):
    stations = glob.glob(Data + '/*')
    stalist = []
    rfcount =[]
    noevents =[]

    for station in stations:
        sta = station.replace(Data+'/MY.','')
        stalist.append(sta)

        file = station + '/selected_RFs_'+noisefilter+filt+'.dat'
        
        
        if os.path.isfile(file):
            with open(file, 'r') as fp:
                x = len(fp.readlines())
                events = glob.glob(station+'/*/*.PICKLE')
        else: 
            x = 0
            events = glob.glob(station+'/*/*/*.PICKLE')
        rfcount.append(x)

        noevents.append(len(events))

    plt.figure(figsize=(12,4))
    plt.bar(np.arange(len(stalist)),rfcount, tick_label=stalist, color='mediumvioletred',alpha=.4, edgecolor='mediumvioletred', label='RFs used')
    plt.bar(np.arange(len(stalist)),noevents, tick_label=stalist, color='coral',alpha=.4, edgecolor='coral', label='Events')
    # plt.yticks(np.arange(0,121, 10))
    plt.ylabel('Frequency')

    plt.legend(frameon=False)

    plt.savefig('selected_RFs_'+noisefilter+filt+'_per_station.png')
    plt.savefig('selected_RFs_'+noisefilter+filt+'_per_station.pdf')