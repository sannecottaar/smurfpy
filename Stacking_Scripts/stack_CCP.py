# common conversion point stacking

# Uses routines from common_converion_point_stack.
# Needs a directory '../CCP_Volumes'


import common_conversion_point_stack as CCP_stack
import glob


#########################
# Set parameters for stack
########################
name= 'CCP_WUS'
rffilter='jgf1'  # RF used
conversion='prem' # conversion use
factor=2. # doubling the fresnel zone width
newstack = True #Starting a new stack (True) or adding data to an old one (False)
lonmin=-130.
lonmax=-100.
lonrez=(lonmax-lonmin)*2.+1. # every 0.5 degrees
latmin=20.
latmax=50.
latrez=(latmax-latmin)*2.+1. # every 0.5 degrees
depmin=60
depmax=1300
deprez=(depmax-depmin)/2. # every 2 km




## Intialize stack or load latest stack
CCP= CCP_stack.ccp_volume()#

if newstack:
    CCP.start_empty_volume(name= name, filter=rffilter, conversion=conversion, factor=factor, lonmin = lonmin, lonmax= lonmax, lonrez=lonrez, latmin=latmin, latmax=latmax, latrez=latrez,depmin=depmin, depmax=depmax, deprez=deprez)
else:
    CCP.load_latest(name= name,filter=rffilter, conversion=conversion, factor=factor)


## Data to add (currently looping through all RFs in station directories, this could be adapted
sta =glob.glob('../Data/*')


for i in range(len(sta)):
    rflist=glob.glob(sta[i]+'/*PICKLE')
    CCP.addlist(rflist,name=name,filter=rffilter, conversion=conversion, factor=factor)


