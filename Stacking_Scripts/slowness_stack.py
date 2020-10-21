'''
Slowness_Stack_Error_New.py
Used to distinguish which signals in the time or depth against epicentral angle stack are converted phases and which signals are multiples.
---Multiples come from shallower angles, therefore the velocity is lower and the slowness is +ve wrt the direct P wave.
---Converted phases come from steeper angles, therefore the velocity is higher and the slowness is -ve wrt the direct P wave.

This script forms a plot of slowness against time, using a specfic epicentral reference distance.

'''

#-------------------------------Set up------------------------------------

#Import required modules
from obspy import read
import matplotlib.pyplot as plt
import os.path
import time
import glob
import numpy as np
import pickle
from shapely.geometry import Point, Polygon

#Array of zeros, slowness stack high and 2*length of RF wide
#201 as 1 to -1 in steps of 0.01
STACK = np.zeros([201,1751])

#Define epicentral distance reference
epi_ref = 60 #Choose based upon epicentral distance range of data, mid-point is good first estimate

#Define epicentral distance bounds
epimin = 30 
epimax = 90 

#Define constraints
depth = 410
lat1 = 89.
lon1 = -179.
lat2 = lat1
lon2 = 179.
lat3 = -89.
lon3 = lon2
lat4 = lat3
lon4 = lon1

box = Polygon([(lat1,lon1),(lat2,lon2),(lat3,lon3),(lat4,lon4)])
#Set some values
filt = 'jgf1'

#Make directory for outputs
savepath='../Slowness_Stacks'
if not os.path.exists(savepath):
    os.makedirs(savepath)
savepath2=savepath+'/Figures/'
if not os.path.exists(savepath2):
    os.makedirs(savepath2)

#-------------------------------Loop through events---------------------------------------

direcs = glob.glob('../Data/*')
stalist = []

#Make list of all events
for direc in direcs:
    print(direc)
    if os.path.isfile(direc + '/selected_RFs_' + filt + '.dat'):
        with open(direc + '/selected_RFs_' + filt + '.dat') as a:
            starfs = a.read().splitlines()
            for line in starfs:
                stalist.append(line)

#Set up checking counts and event lists
count_yes = 0
count_no = 0
c = 0

#Loop through events
for i in range(len(stalist)):


    print (stalist[i])
    c = c +1
    print (c)

    #Read in the file
    seis = read(stalist[i], format = 'PICKLE')
    RF = getattr(seis[0],filt)

    #--------Constrain for epicentral distance and pierce point lat/lon-----------

    #Extract epicentral distance of trace
    epi_RF = seis[0].stats['dist']
    stel = round(seis[0].stats['stel'],3) #  Use above ground staitons for simplicity.
    #Make RFs the same length
    while len(RF['iterativedeconvolution']) < 1751:
        RF['iterativedeconvolution'] = np.append(RF['iterativedeconvolution'],0)
        RF['time'] = np.append(RF['time'],RF['time'][-1]+0.1)

    #Extract pierce point lat/lon for Pds at depth, d
    trace_pierce = seis[0].stats['piercepoints']['P'+str(depth)+'s'][str(depth)]

    #Define lat and lon
    latpp = float(trace_pierce[1])
    lonpp = float(trace_pierce[2])

    #Make the few positive (east) lon degrees an extension of negative
    if lonpp >= 0:
        lonpp = lonpp - 360

    #Make Shapely points
    point = Point(float(trace_pierce[1]), float(trace_pierce[2]))

    #Check whether event has epicentral distance and piercepoint lat,lon needed
    if epimin <= epi_RF <= epimax and stel >= 0.0:
        if box.contains(point):
            print ('yes')
            count_yes = count_yes + 1


            #--------------------Make an array of slowness----------------

            #List values between -1.00 and +1.00 in steps of 0.01 - slowness in s/deg

            #Make a vector of integers for the slowness
            slow_int = range(-100,101,1)

            #Divide the integers to get the exact slowness values
            slow = [x / 100. for x in slow_int]
            #print(len(slow))

            #Loop through each slowness value
            for j in slow_int:
                s = slow[j]

                #Shift the epicentral distance in relation to the reference
                epi_dist = epi_RF - epi_ref

                #Calculate delta t (timeshift)
                timeshift = s * epi_dist


                #--------Convert the RF into the right format--------------

                #Extract the amplitude of the RF
                RF_amp = getattr(seis[0],filt)['iterativedeconvolution']

                #Create line of zeros length of 40s on either side of the RF (so shifted traces can be stacked)
                #40s either side is  80/float(seis[0].stats['delta']), so 80/sampling time (0.1s)
                RFtemp = np.zeros(len(getattr(seis[0],filt)['iterativedeconvolution'])+int(80./float(seis[0].stats['delta'])))

                #Take the middle section of 0s of RFtemp and make it equal to RF_amp
                RFtemp[int(40/float(seis[0].stats['delta'])):int(40/float(seis[0].stats['delta']))+len(getattr(seis[0],filt)['iterativedeconvolution'])] = RF_amp


                #--------Set up slowness stack-------------------------

                #Set up timeshift of sample
                timeshift_sample = int(timeshift/float(seis[0].stats['delta']))

                #Add the timeshifted RF to the stack
                STACK[j,:] = STACK[j,:] + RFtemp[int(40/float(seis[0].stats['delta']) + timeshift_sample) : int(40/float(seis[0].stats['delta'])) + timeshift_sample + len(getattr(seis[0],filt)['iterativedeconvolution'])]


        #Add to count if not in bounds for lat/lon
        else:
            print ('no')
            count_no = count_no + 1

    #Add to count if not in bounds for epicentral distance
    else:
        print ('no')
        count_no = count_no + 1


#------------------Standard Error Calculation--------------------------

#Reset Counter
c = 0

#Find average of STACK
STACK = STACK / count_yes

#Set up blank error stack
STACK_ERR = np.zeros([201,1751])

#Loop through events
for i in range(len(stalist)):
    seis = read(stalist[i], format = 'PICKLE')
    RF = getattr(seis[0],filt)
    print (stalist[i])
    c = c +1
    print (c)

    #--------Constrain for epicentral distance and pierce point lat/lon-----------

    #Make RFs the same length
    while len(RF['iterativedeconvolution']) < 1751:
        RF['iterativedeconvolution'] = np.append(RF['iterativedeconvolution'],0)
        RF['time'] = np.append(RF['time'],RF['time'][-1]+0.1)

    #Extract epicentral distance of trace
    epi_RF = seis[0].stats['dist']
    stel = round(seis[0].stats['stel'],3) #  Use above ground staitons for simplicity.

    #Extract pierce point lat/lon for Pds at depth, d
    trace_pierce = seis[0].stats['piercepoints']['P'+str(depth)+'s'][str(depth)]

    #Define lat and lon
    latpp = float(trace_pierce[1])
    lonpp = float(trace_pierce[2])

    #Make the few positive (east) lon degrees an extension of negative
    if lonpp >= 0:
        lonpp = lonpp - 360

    #Make Shapely points
    point = Point(float(trace_pierce[1]), float(trace_pierce[2]))

    #Check whether event has epicentral distance and piercepoint lat,lon needed
    if epimin <= epi_RF <= epimax and stel >= 0.0:
        if box.contains(point):
            print ('yes')

            #--------------------Make an array of slowness----------------

            #List values between -1.00 and +1.00 in steps of 0.01 - slowness in s/deg

            #Make a vector of integers for the slowness
            slow_int = range(-100,101,1)

            #Divide the integers to get the exact slowness values
            slow = [x / 100. for x in slow_int]

            #Loop through each slowness value
            for j in slow_int:
                s = slow[j]

                #Shift the epicentral distance in relation to the reference
                epi_dist = epi_RF - epi_ref

                #Calculate delta t (timeshift)
                timeshift = s * epi_dist


                #--------Convert the RF into the right format--------------

                #Extract the amplitude of the RF
                RF_amp = getattr(seis[0],filt)['iterativedeconvolution']

                #Create line of zeros length of 40s on either side of the RF (so shifted traces can be stacked)
                #40s either side is  80/float(seis[0].stats['delta']), so 80/sampling time (0.1s)
                RFtemp = np.zeros(len(getattr(seis[0],filt)['iterativedeconvolution'])+int(80/float(seis[0].stats['delta'])))

                #Take the middle section of 0s of RFtemp and make it equal to RF_amp
                RFtemp[int(40/float(seis[0].stats['delta'])):int(40/float(seis[0].stats['delta']))+len(getattr(seis[0],filt)['iterativedeconvolution'])] = RF_amp


                #--------Set up slowness stack-------------------------

                #Set up timeshift of sample
                timeshift_sample = int(timeshift/float(seis[0].stats['delta']))

                #Find timeshifted RF
                RF_timeshift = RFtemp[int(40/float(seis[0].stats['delta'])) + timeshift_sample : int(40/float(seis[0].stats['delta'])) + timeshift_sample + len(getattr(seis[0],filt)['iterativedeconvolution'])]

                #Find deviation from the mean and add to error stack
                STACK_ERR[j,:] = STACK_ERR[j,:] + ((STACK[j,:] - RF_timeshift)**2)

print(count_yes)

#Find the standard deviation by dividing by the sample number and square root
SD = np.sqrt(STACK_ERR/(count_yes-1))

#Find the standard error by dividing the SD by the square root of the sample size
SE = SD/(np.sqrt(count_yes))

#Extract time axis
time=getattr(seis[0],filt)['time']



#-----------------Create Stack only within errors-------------

#Create copy of STACK
STACK2 = STACK.copy()

#Loop
for c in range(201):
    for d in range(1751):
        if STACK2[c,d] > 0:
            STACK2[c,d] = STACK2[c,d] - 2*SE[c,d]
            if STACK2[c,d] < 0:
                STACK2[c,d] = 0
        elif STACK2[c,d] < 0:
            STACK2[c,d] = STACK2[c,d] + 2*SE[c,d]
            if STACK2[c,d] > 0:
                STACK2[c,d] = 0



#-------------Do the same for this plot--------------

#Normalisation value
NORMALIZATION =  1.2

#Find highest value of the STACK and make this 1
STACK2 = STACK2 / np.max(np.abs(STACK2))

#This will bring out the smaller Pds phases that have a smaller amplitude
STACK2 = STACK2/NORMALIZATION

#Print checks
print (count_yes)
print (count_no)

#Set size of image
fig = plt.figure(figsize=(12,4))


#---------------------Plot Predicted Travel Times and Slowness of converted phases-----------------

#Read in the predicted travel time and slowness differences from the appropriate file
data = np.genfromtxt('../Tools/Travel_Times_Slowness/TTS_Pds_'+str(epi_ref)+'.dat', delimiter='\t')

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[:,1]
y = data[:,2]

#Plot the curves on the stack
plt.plot(x,y)


#---------------------Plot Predicted Travel Times and Slowness of PPvdp-----------------

#Read in the predicted travel time and slowness differences from the appropriate file
data = np.genfromtxt('../Tools/Travel_Times_Slowness/TTS_PPvdp_'+str(epi_ref)+'.dat', delimiter='\t')

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[:,1]
y = data[:,2]

#Plot the curves on the stack
plt.plot(x,y)
#plt.annotate('PPvdp',xy=(130,0.63))


#---------------------Plot Predicted Travel Times and Slowness of PPvds-----------------

#Read in the predicted travel time and slowness differences from the appropriate file
data = np.genfromtxt('../Tools/Travel_Times_Slowness/TTS_PPvds_'+str(epi_ref)+'.dat', delimiter='\t')

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[:,1]
y = data[:,2]

#Plot the curves on the stack
plt.plot(x,y)
#plt.annotate('PPvds',xy=(137,0.17))


#---------------------Plot Predicted points for discontinuities-----------------

#Read in the predicted travel time and slowness differences from the appropriate file
data = np.genfromtxt('../Tools/Travel_Times_Slowness/TTS_Pds_'+str(epi_ref)+'.dat', delimiter='\t')

#220km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[1,1]
y = data[1,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#310km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[2,1]
y = data[2,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#410km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[3,1]
y = data[3,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#550km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[4,1]
y = data[4,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#660km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[5,1]
y = data[5,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#971km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[6,1]
y = data[6,2]

#1000km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[7,1]
y = data[7,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#1071km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = data[8,1]
y = data[8,2]

#-----------------------PPvdp

#Read in the predicted travel time and slowness differences from the appropriate file
datap = np.genfromtxt('../Tools/Travel_Times_Slowness/TTS_PPvdp_'+str(epi_ref)+'.dat', delimiter='\t')
#220km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datap[1,1]
y = datap[1,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#310km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datap[2,1]
y = datap[2,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#410km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datap[3,1]
y = datap[3,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#550km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datap[4,1]
y = datap[4,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#660km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datap[5,1]
y = datap[5,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#-----------------------PPvds

#Read in the predicted travel time and slowness differences from the appropriate file
datas = np.genfromtxt('../Tools/Travel_Times_Slowness/TTS_PPvds_'+str(epi_ref)+'.dat', delimiter='\t')
#220km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datas[1,1]
y = datas[1,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#310km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datas[2,1]
y = datas[2,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#410km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datas[3,1]
y = datas[3,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#550km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datas[4,1]
y = datas[4,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#660km

#Get the travel time difference, set as x, and the slowness difference, set as y
x = datas[5,1]
y = datas[5,2]

#Plot the curves on the stack
plt.plot(x,y, marker="x", color='black', markersize=8, markeredgewidth=1.6)

#-----------------Plot the stack----------------------------

#Plotting using contourf (time_axis (RF time), slowness_range, slowness_matrix, contour_levels , colourmap , extend (extend grid if data lies outside range))
contour_levels = np.arange(-.3,.3,.02)  #(-1.,1.01,.1) (min,max,interval)
cs = plt.contourf(time,slow,STACK2,contour_levels, cmap=plt.cm.seismic, extend="both")

#Set axes limits
plt.ylim(-0.8,0.8)
plt.xlim(0,150)

#Axes labels and title
plt.ylabel ('Slowness (s/deg)')
plt.xlabel ('Time (s)')
plt.suptitle('Slowness Stack - Epicentral Ref. '+str(epi_ref)+' - '+str(depth)+'km - Filt '+str(filt)+' \n No. of RFs: '+str(count_yes)+' Lat/Lon '+'North '+'  Norm: '+str(NORMALIZATION))

#Save figures
plt.savefig(savepath+'/Figures/Slowness_'+str(filt)+'_'+str(count_yes)+'RFs'+'.pdf')
plt.savefig(savepath+'/Figures/Slowness_'+str(filt)+'_'+str(count_yes)+'RFs'+'.png')

#Show the plot
plt.show()

#-------------------------------Save stack as pickle file----------------------------

outpickle=dict()
outpickle['STACK']=STACK2
outpickle['STACK_nonnorm']=STACK
outpickle['STACK_SE']=SE
outpickle['time']=time
outpickle['slow']=slow
outpickle['PP_depth']=depth
outpickle['no_RFs']=count_yes
outpickle['latmin']=lat1
outpickle['latmax']=lat3
outpickle['lonmin']=lon1
outpickle['lonmax']=lon2
outpickle['filter']=filt
outpickle['NORMALIZATION']=NORMALIZATION

output_file=str(savepath+'/Slowness_'+str(filt)+'_'+str(count_yes)+'RFs'+'.PICKLE')
out_put_write=open(output_file,'wb')
pickle.dump(outpickle, out_put_write)
out_put_write.close
