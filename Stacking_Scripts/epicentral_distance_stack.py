#
import obspy
from obspy.core import trace
from obspy import read

import matplotlib.pyplot as plt
import glob
import numpy as np
import pickle

comp = 'jgf1'

stadirs = glob.glob('../Data/*')

start = 30.
stop = 90.
n = 13
step = (stop - start)/(n-1)
start_bins = start + ((stop - start) / n)

deg = np.linspace(start, stop, n)
#print(deg)
degplot = []

for i in range(len(deg)):
    degplot.append(deg[i] + (step / 2))

#print(degplot)

stack = np.zeros((len(deg), 1751))
stacked = np.transpose(stack)

count = [[] for j in range(n)]
for j in range(n):
    count[j] = 0




for stadir in stadirs:

    # loop through events
    evntlist=glob.glob(stadir+'/*.PICKLE')

    # Loop through data
    if(len(evntlist)>0):
        for i in range(len(evntlist)): #range(cat.count()):
                seis=read(evntlist[i],format='PICKLE')
                distdg=seis[0].stats['dist']
                #print(distdg)
                RF=trace.Trace()
                #print(RF)

                # Select epicentral distance bin
                for bin in range(len(deg)):
                    if distdg <= degplot[bin]:
                        k = bin

                        break

                # Get Rf data
                if comp=='rff2':
                    RF.data= np.real(seis[0].rff2['iterativedeconvolution'])
                if comp=='rff1':
                    RF.data= np.real(seis[0].rff1['iterativedeconvolution'])
                if comp=='jgf1':
                    RF.data= np.real(seis[0].jgf1['iterativedeconvolution'])
                if comp=='jgf2':
                    RF.data= np.real(seis[0].jgf2['iterativedeconvolution'])

                # Stack data
                stack[k,] += RF.data
                count[k] +=1

# Take mean
for j in range(n):
    stack[j,] = stack[j,] / count[j]

transpose = np.transpose(stack)


# Plot data
contour_levels = np.arange(-.5,.51,.025)
plt.contourf(transpose,contour_levels,rasterized=True, cmap='seismic',extend='both')
plt.xlabel('Epicentral distance (deg)')
plt.ylabel('Time (s)')
plt.xticks(np.arange(n), np.linspace(start, stop, num=n, endpoint=True))
plt.yticks(np.linspace(0, 1751, num=8), np.linspace(-25, 150, num=8, endpoint=True))
plt.title('Epicentral distance stack')

plt.show()

for i in range(n):
    for j in range(len(RF)):
        stacked[j,i] = transpose[j,i] + i

plt.plot(stacked)
plt.ylabel('Epicentral distance (deg)')
plt.xlabel('Time (s)')
plt.title('Receiver functions at each epicentral distance')
plt.yticks(np.arange(n), np.linspace(start, stop, num=n, endpoint=True))
plt.xticks(np.linspace(0, 1751, num=8), np.linspace(-25, 150, num=8, endpoint=True))

plt.show()
