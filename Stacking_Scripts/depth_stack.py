'''
depth_stack.py
written by Matthew Kemp and Sanne Cottaar

Stacks all the RFs within the bounds for the depth stated producing one trace.
Normalised at various points to emphasise the smaller phases.
Reruns through all data to determine standard error.
Adaptable to plot the depths of various conversions.

To run this stack, conversions need to be computed first. 
For PREM conversions use Migration_Scripts/convert_to_depth_obspy.py


'''
#--------------------------Set Up----------------------

#Import all the relevant modules
import obspy
from obspy import read
from obspy.core import Stream
from obspy.core import trace
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os.path
import time
import glob
import shutil
import numpy as np
from obspy import UTCDateTime
import subprocess
import pickle
import sys
from scipy import interpolate


#Set figure parameters
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (8,10),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

#Create depth stack
depth_space = np.linspace(0,1200,1201)
STACK=np.zeros(len(depth_space))

#Set some values
direc = '../Data'
filt = 'jgf1' # RF type to use
conversion = 'prem' # Converstion to use

norm_depth = 150 # depth beyond which to normalize
norm_fact = 0.2 # Normalization factor

savepath='./'

#Define constraints
dp = 410
latmin = -25
latmax = 0
lonmin = -178
lonmax = -160

#Make list of all events
stalist = glob.glob(direc+'*/*/*.PICKLE')

#Set up checking counts and event lists
count_yes = 0
count_no = 0
List_yes = []
List_no = []
c = 0

#--------------------------Main Loop----------------------

#Loop through events
for s in range(len(stalist)):   
    print (stalist[s])
    c = c +1
    print (c)

    #Read in the file
    seis = read(stalist[s],format='PICKLE')
    RF = getattr(seis[0],filt)
    
    #Extract the pierce points
    trace_pierce = seis[0].stats['piercepoints']['P'+str(dp)+'s'][str(dp)]
    print (trace_pierce)

    #Define lat and lon of event pierce points
    latpp = float(trace_pierce[1])
    lonpp = float(trace_pierce[2])
    
    #Make the few positive (east) lon degrees an extension of negative
    if lonpp >= 0:
        lonpp = lonpp - 360

    #Check whether event pierce points located within bounds of regions
    if (latmin <= latpp <= latmax) and (lonmin <= lonpp <= lonmax):
        print ('yes')

        #Add to list and count
        List_yes.append(stalist[s])
        count_yes = count_yes + 1


    
        # for larger data sets, we should consider writing out the interpolated values...
        interp = interpolate.interp1d(seis[0].conversions[conversion]['depthsfortime'],RF['iterativedeconvolution'] )
        RFdata_depth = interp(depth_space)
            


        #Add this amplitude to the stack
        STACK[:] = STACK[:] + RFdata_depth[:]

    #If not within bounds of region
    else:
        print ('no')

        #Add to list and count
        List_no.append(stalist[s])
        count_no = count_no + 1

#Extract the regular depth intervals from the event
depth = depth_space

#Find average amplitudes
amp_rel = STACK[:]/count_yes


#--------------------------Error calculation----------------------

#Set up total for deviations
dev_sq_tot = 0
c = 0

#Loop through events
for s in range(len(stalist)):
    print (stalist[s])
    c = c +1
    print (c)
    seis = read(stalist[s],format='PICKLE')
    trace_pierce = seis[0].stats['piercepoints']['P'+str(dp)+'s'][str(dp)]
    print (trace_pierce)
    latpp = float(trace_pierce[1])
    lonpp = float(trace_pierce[2])
    if lonpp >= 0:
        lonpp = lonpp - 360
    if (latmin <= latpp <= latmax) and (lonmin <= lonpp <= lonmax):
        print ('yes')

        interp = interpolate.interp1d(seis[0].conversions[conversion]['depthsfortime'],RF['iterativedeconvolution'])
        RFdata_depth = interp(depth_space)

        amp = RFdata_depth
 
	
        #Find deviation from the mean
        dev = amp_rel[:] - amp[:]
        
        #Sqaure the deviation to make absolute
        dev_sq = dev**2
	
        #Add this to total
        dev_sq_tot = dev_sq_tot + dev_sq

    #If not within bounds of region
    else:
        print ('no')

#Find the standard deviation by dividing by the sample numberand square root
SD = np.sqrt(dev_sq_tot/(count_yes-1))

#Find the standard error by dividing the SD by the square root of the sample size
SE = SD/(np.sqrt(count_yes))

#Normalise max value (direct P) to 1
amp_rel[:]=amp_rel[:]/(np.nanmax(np.abs(amp_rel[:])))

#Attempts to find peaks automatically
#import scipy
#from scipy import signal
#scipy.signal.find_peaks_cwt(amp_rel, np.arange(1,100))

#Add and subtract SE from amp_rel to show spread of error
amp_rel_1SE_P = amp_rel[:] + SE
amp_rel_1SE_N = amp_rel[:] - SE
amp_rel_2SE_P = amp_rel[:] + SE*2
amp_rel_2SE_N = amp_rel[:] - SE*2


int_depth = np.argmin(np.abs(depth_space-norm_depth))
#Re-normalise later sections to emphasise the small phases
amp_rel[int_depth:-1]=amp_rel[int_depth:-1]/norm_fact


#Re-normalise later sections for errors
amp_rel_1SE_P[int_depth:-1]=amp_rel_1SE_P[int_depth:-1]/norm_fact
amp_rel_1SE_N[int_depth:-1]=amp_rel_1SE_N[int_depth:-1]/norm_fact
amp_rel_2SE_P[int_depth:-1]=amp_rel_2SE_P[int_depth:-1]/norm_fact
amp_rel_2SE_N[int_depth:-1]=amp_rel_2SE_N[int_depth:-1]/norm_fact


#Print checks
#print (List_yes)
#print (List_no)
print (count_yes)
print (count_no)


#--------------------------Plotting----------------------

#Plot the relative amplitudes against the depth points
plt.plot(amp_rel, depth, 'k', lw = 1.5)

#Plot error lines
plt.plot(amp_rel_1SE_P, depth, 'k--', lw = 0.75, color='0.75')
plt.plot(amp_rel_1SE_N, depth, 'k--', lw = 0.75, color='0.75')
plt.plot(amp_rel_2SE_P, depth, 'k--', lw = 1, color='0.5')
plt.plot(amp_rel_2SE_N, depth, 'k--', lw = 1, color='0.5')

#Fill the postive values with red, negative with blue
plt.fill_betweenx(depth, 0, amp_rel_2SE_N, where=amp_rel_2SE_N >= 0, facecolor=[1, 0.4, 0.4])
plt.fill_betweenx(depth, 0, amp_rel_2SE_P, where=amp_rel_2SE_P <= 0, facecolor=[0.4, 0.4, 1])

#Make label for y axis (none for x)
plt.ylabel('depth (km)', fontsize=16)
plt.xticks([])

#Create dotted lines to show where there are new normalisations
x=[-0.5,1.1]
y1=[150,150]
plt.plot(x,y1, 'k--')

#Set axis limits and invert y axis
plt.gca().set_xlim([-0.5,1.1])
plt.gca().set_ylim([0,1200])
plt.gca().invert_yaxis()
plt.gca().tick_params(labelsize=16)

#Make title of plot
plt.suptitle('Depth Stack - PP: '+str(dp)+'km - No. of RFs: '+str(count_yes)+'\n Lat/Lon: '+str(latmin)+'/'+str(latmax)+'/'+str(lonmin)+'/'+str(lonmax)+'\n Filter: '+str(filt)+' Conversion: '+str(conversion))

#Save figures
plt.savefig(savepath +'/Depth_Stack_'+str(filt)+'_'+str(conversion)+'_'+str(count_yes)+'RFs_loc'+str(latmin)+'_'+str(latmax)+'_'+str(lonmin)+'_'+str(lonmax)+'.pdf')
plt.savefig(savepath+'/Depth_Stack_'+str(filt)+'_'+str(conversion)+'_'+str(count_yes)+'RFs_loc'+str(latmin)+'_'+str(latmax)+'_'+str(lonmin)+'_'+str(lonmax)+'.png')

#Show the plot
plt.show()

#-------------------------------Save stack as pickle file----------------------------

#epi_range,time,STACK_NEW,counter

outpickle=dict()
outpickle['depth']=depth
outpickle['amp_rel']=amp_rel
outpickle['amp_rel_1SE_P']=amp_rel_1SE_P
outpickle['amp_rel_1SE_N']=amp_rel_1SE_N
outpickle['amp_rel_2SE_P']=amp_rel_2SE_P
outpickle['amp_rel_2SE_N']=amp_rel_2SE_N
outpickle['PP_depth']=dp
outpickle['no_RFs']=count_yes
outpickle['latmin']=latmin
outpickle['latmax']=latmax
outpickle['lonmin']=lonmin
outpickle['lonmax']=lonmax
outpickle['filter']=filt
outpickle['conversion']=conversion

output_file=str(savepath+'/Depth_Stack_'+str(filt)+'_'+str(conversion)+'_'+str(count_yes)+'RFs_loc'+str(latmin)+'_'+str(latmax)+'_'+str(lonmin)+'_'+str(lonmax)+'.PICKLE')
out_put_write=open(output_file,'wb')
pickle.dump(outpickle, out_put_write)
out_put_write.close
