
from obspy import read
import matplotlib.pyplot as plt
import glob
import sys

# Command line help
if len(sys.argv) != 2 or str(sys.argv[1]).lower() == 'help':
    print('\n')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print(sys.argv[0])
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('Description:           [OPTIONAL] Plots the perstation distribution of "Acceptable - Green" and "Removed - red" events')
    print('                       as a funciton of EQ magnitude and epicentral distance.')
    print('Inputs:                Data directory (usually ../Data/)')
    print('Outputs:               On-screen plotting')
    print('Usage:                 >> python3  7_plot_data_selection.py filterband')
    print('Options [1]:           jgf1, jgf2, jgf3, tff1, tff2, tff3, tff4 or tff5')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('\n')
    sys.exit()

filt=str(sys.argv[1])

stadirs = glob.glob('../Data/*')

mags = []
dist =[]
for stadir in stadirs:
    print(stadir)

    # loop through events
    stalist=glob.glob(stadir+'/*.PICKLE')

    station_file = open(stadir + '/selected_RFs_'+str(filt)+'.dat','r')
    stalist=station_file.read().strip().split()
    print(stalist)
    # Loop through data
    if(len(stalist)>0):
        for i in range(len(stalist)): #range(cat.count()):
                print(stalist[i])
                seis=read(stalist[i],format='PICKLE')
                mags.append(seis[0].stats['event']['magnitudes'][0]['mag'])
                dist.append(seis[0].stats['dist'])

plt.scatter(mags,dist,s=15,c='g', edgecolor='none')


mags = []
dist =[]
for stadir in stadirs:
    print(stadir)

    # loop through events
    stalist_all=glob.glob(stadir+'/*.PICKLE')
    station_file = open(stadir + '/selected_RFs_'+str(filt)+'.dat','r')
    stalist_good=station_file.read().strip().split()
    stalist_poor = [x for x in stalist_all if x not in stalist_good]
    print(stalist_poor)
    # Loop through data
    if(len(stalist_poor)>0):
        for i in range(len(stalist_poor)): #range(cat.count()):
                seis=read(stalist_poor[i],format='PICKLE')
                mags.append(seis[0].stats['event']['magnitudes'][0]['mag'])
                dist.append(seis[0].stats['dist'])

plt.scatter(mags,dist,s=5,c='r', edgecolor='none')

plt.show()

                
 


                
 
