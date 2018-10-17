#

from obspy import read
import matplotlib.pyplot as plt
import time
import glob
import numpy as np
from obspy import UTCDateTime

direc = 'DataRF'
flag = 'SV'
filt = 'jgf1'



stadirs = glob.glob(direc+'/*')


for stadir in stadirs:
    print(stadir)
    with open(stadir+'/selected_RFs_jgf1.dat','r') as f:
        goodrfs= f.read().replace('\n', '')

    # loop through events
    stalist=glob.glob(stadir+'/*.PICKLE')
    print(stalist)
    c=0
    # Loop through data
    if(len(stalist)>0):
        for i in range(len(stalist)): #range(cat.count()):
                print(stalist[i])
                seis=read(stalist[i],format='PICKLE')
 
                distdg=seis[0].stats['dist']
                if stalist[i] in goodrfs:
                    good=True
                    print('YAY',seis[0].stats['event'].magnitudes[0].mag)
                else:
                    good=False
                    print('NO',seis[0].stats['event'].magnitudes[0].mag)


                tshift=UTCDateTime(seis[0].stats['starttime'])-seis[0].stats['event'].origins[0].time
  
                #Ptime=Ptime
                plt.subplot(1,3,1)
                vertical = seis.select(channel='*HZ')[0]
                vertical.filter('bandpass', freqmin=0.01,freqmax=.1, corners=2, zerophase=True)
                windowed=vertical[np.where(vertical.times()>seis[0].stats.traveltimes['P']-100) and np.where(vertical.times()<seis[0].stats.traveltimes['P']+100)]
                norm=np.max(np.abs(windowed))
                if good:
                    plt.plot(vertical.times()-seis[0].stats.traveltimes['P'], vertical.data/norm+np.round(distdg),'k')
                else:
                    plt.plot(vertical.times()-seis[0].stats.traveltimes['P'], vertical.data/norm+np.round(distdg),'r')                    
                #plt.plot(seis[0].stats.traveltimes['P'],np.round(distdg),'.b')
                #plt.plot(seis[0].stats.traveltimes['S'],np.round(distdg),'.g')
                plt.xlim([-25,150])
                plt.ylim([30,92])
                plt.subplot(1,3,2)
                radial = seis.select(channel='*HR')[0]

                radial.filter('bandpass', freqmin=0.01,freqmax=.1, corners=2, zerophase=True)
                windowed=vertical[np.where(radial.times()>seis[0].stats.traveltimes['P']-100) and np.where(radial.times()<seis[0].stats.traveltimes['P']+100)]
                norm=np.max(np.abs(windowed))
                if good:
                    plt.plot(radial.times()-seis[0].stats.traveltimes['P'], radial.data/norm+np.round(distdg),'k')
                else:
                    plt.plot(radial.times()-seis[0].stats.traveltimes['P'], radial.data/norm+np.round(distdg),'r')
                plt.xlim([-25,150])                  
                plt.plot(seis[0].stats.traveltimes['P'],np.round(distdg),'.b')
                plt.plot(seis[0].stats.traveltimes['S'],np.round(distdg),'.g')
                plt.ylim([30,92])
            
                plt.subplot(1,3,3)
                RF=getattr(seis[0],filt)['iterativedeconvolution']
                time=getattr(seis[0],filt)['time']
                if good:
                    plt.plot(time, RF/np.max(np.abs(RF))+np.round(distdg),'k')
                else:
                    plt.plot(time, RF/np.max(np.abs(RF))+np.round(distdg),'r')                   

                
                
                    
        plt.subplot(1,3,1)
        plt.title('vertical')
        plt.ylabel('distance')
        plt.xlabel('time')
        plt.subplot(1,3,2)
        plt.title('radial')
        plt.ylabel('distance')
        plt.xlabel('time')
        plt.subplot(1,3,3)
        plt.title('receiver functions')
        plt.ylabel('distance')
        plt.xlabel('time')
        #plt.xlim([-150,1000])        
        plt.show()






