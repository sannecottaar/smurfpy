import matplotlib.pyplot as plt
import numpy as np
import glob
import os

#Command line help
if len(sys.argv) != 2 or str(sys.argv[1]).lower() == 'help':
    print('\n')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print(sys.argv[0])
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('Description:           [OPTIONAL] Plot how many RFs are generated and selected per station as a bar chart.')
    print('Inputs:                Data directory (usually ../Data/), filter band')
    print('Outputs:               n/a\n')
    print('Usage:                 >> python3 RFsPerStation.py datadirectory filterband')
    print('Options [1]:           jgf1, jgf2, jgf3, tff1, tff2, tff3, tff4 or tff5 [str]')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('\n')
    sys.exit()

def rfsperstation(Data, filt):
    stations = glob.glob(Data + '/*')
    stalist = []
    rfcount =[]
    noevents =[]

    for station in stations:
        sta = station.replace(Data+'/MY.','')
        stalist.append(sta)

        file = station + '/selected_RFs_'+filt+'.dat'
        
        
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

    plt.savefig('selected_RFs_'+filt+'_per_station.png')
    plt.savefig('selected_RFs_'+filt+'_per_station.pdf')
    
Data = sys.argv[1]
filt = sys.argv[2]
rfsperstation(Data, filt)
