# common conversion point stacking

# Uses routines from common_converion_point_stack.
# Needs a directory '../CCP_Volumes'


import common_conversion_point_stack as CCP_stack
import glob
import os
import time
import sys
from obspy import read

#########################
# Set parameters for stack
########################
name= 'CCP_South_Africa'
rffilter='jgf1'  # RF used
conversion='prem' # conversion use
factor=2. # doubling the fresnel zone width
newstack = True #Starting a new stack (True) or adding data to an old one (False)
lonmin=11.
lonmax=40.
latmin=-36.
latmax=-14.
depmin=60.
depmax=1300.
lonrez=(lonmax-lonmin)*2.+1. # every 0.5 degrees
latrez=(latmax-latmin)*2.+1. # every 0.5 degrees
deprez=(depmax-depmin)/2. # every 2 km

script_start = time.time()

## Intialize stack or load latest stack
CCP= CCP_stack.ccp_volume()#

if newstack:
    CCP.start_empty_volume(name= name, filter=rffilter, conversion=conversion, factor=factor, lonmin = lonmin, lonmax= lonmax, lonrez=lonrez, latmin=latmin, latmax=latmax, latrez=latrez,depmin=depmin, depmax=depmax, deprez=deprez)
else:
    CCP.load_latest(name= name,filter=rffilter, conversion=conversion, factor=factor)


## Data to add (currently looping through all RFs in station directories, this could be adapted
stations=glob.glob('../Data/*')


for sta in stations:
    rflist=[]
    if os.path.exists(sta + '/selected_RFs_'+str(rffilter)+'.dat'):
        station_file = open(sta + '/selected_RFs_'+str(rffilter)+'.dat','r')
        file_list=station_file.read().strip().split()
        for file in file_list:
            seis=read(file, format='PICKLE')
            if hasattr(seis[0],rffilter):
                rflist.append(file)
        CCP.addlist(rflist,name=name,filter=rffilter, conversion=conversion, factor=factor)

        print('Finished adding RFs... now cleaning...')
        line=open('../CCP_volumes/'+name+'_'+rffilter+'_'+conversion+'_'+str(factor)+'/filenames.dat','r').readlines()
        if len(line)>1:
            rm_file=line[-2].split()[1]
            try:
                print('deleting '+str(rm_file)+'...')
                os.remove(rm_file)
            except:
                print('Failed to delete '+str(rm_file))


print('CCP stacking for '+str(name)+'_'+str(rffilter)+'_'+str(conversion)+'_'+str(factor)+' is COMPLETE!!!')
print('It took', time.time()-script_start, 'seconds.')
