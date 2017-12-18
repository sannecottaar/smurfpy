
from obspy import read
import matplotlib.pyplot as plt
import time
import glob
import numpy as np
from obspy import UTCDateTime
import os



stadirs = glob.glob('../Data/*')

mags = []
dist =[]
for stadir in stadirs:
    print(stadir)

    # loop through events
    stalist=glob.glob(stadir+'/*.PICKLE')
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
    stalist=glob.glob(stadir+'/auto_removed/*.PICKLE')
    print(stalist)
    # Loop through data
    if(len(stalist)>0):
        for i in range(len(stalist)): #range(cat.count()):
                seis=read(stalist[i],format='PICKLE')
                mags.append(seis[0].stats['event']['magnitudes'][0]['mag'])
                dist.append(seis[0].stats['dist'])

plt.scatter(mags,dist,s=5,c='r', edgecolor='none')

plt.show()

                
 


                
 
