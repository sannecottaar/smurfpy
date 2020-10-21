'''
epicentral_distance_stack.py
Scripted by Matthew Kemp, Annemijn van Stiphout, Kieran Gilmore, others?

Stacks the RFs in bins of epicentral distance to show the most prominent features.
Plots a graph of epicentral distance against time, with a red-blue colour map showing positive and negative amplitudes.
The script will ask for various inputs:

  step (width of epicentral distance bin) i.e. 2
  smoothing (adds adjacent bins to smooth the features)
  histogram (plots epicentral distance histogram; how many RFs in each bin, pre-smoothing)
  depth (depth of pierce points and Pds)
  latlon (give min and max for lat and lon constaints on piercepoints)
  predicted travel time (plots the travel time curves for all phases written here) i.e. P1000s

'''


#----------------------------------------Import all the relevant modules--


from obspy import read
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import glob,os
import numpy as np
from shapely.geometry import Point, Polygon
from matplotlib import gridspec

file_char_len = -51
#----------------------------------------User Interface-------------------

# Initial options

bin_size = 5      # in degrees
smoothing = False # adds in neighbouring bins
plot_histogram = True  # plot histogram of data covarage

# Get depth of piercepoints and Pds phase to select data by

depth = 410

# Lat and Lon constraints, piercepoints within this box will be accepted 
lat1 = 89.
lon1 = -179.
lat2 = lat1
lon2 = 179.
lat3 = -89.
lon3 = lon2
lat4 = lat3
lon4 = lon1

box = Polygon([(lat1,lon1),(lat2,lon2),(lat3,lon3),(lat4,lon4)])

# Plotting Travel Time Curves
plot_travel_time_curves = True
phases_to_plot  = ['PP', 'P410s', 'P660s']


#----------------------------------------Set up the arrays for the STACK a

# Define min/max epi dist and no. steps
step = bin_size
epi_steps = (60//step)+1
min_epi = 30
max_epi = 90
epi_range = np.linspace(min_epi, max_epi, epi_steps)

# Set up stack ( matrix [len of RF x no of bins ])
STACK = np.zeros([1751, epi_steps-1])

# Set up counter to see how many RFs are added to each bin
counter = np.zeros([epi_steps-1])

savepath='../Epicentral_Distance_Stack/'
if not os.path.exists(savepath):
    os.makedirs(savepath)

#----------------------------------------Loop through the RF data---------

# Set up checking counts
count_yes = 0
count_no = 0
countwrong = 0

# Set up
direc = '../Data'
filt = 'jgf1'
stadirs = glob.glob(direc + '/*')

# Loop through stations
for stadir in stadirs:

    stalist = []
    if os.path.isfile(stadir + '/selected_RFs_' + filt + '.dat'):
        with open(stadir + '/selected_RFs_' + filt + '.dat') as a:
            starfs = a.read().splitlines()
            for line in starfs:
                stalist.append(line)

    # Loop through events
    for i in range(len(stalist)):
        print (stalist[i])

        # Read in RF
        seis = read(stalist[i], format='PICKLE')

        # Define part RF part of event
        out = getattr(seis[0], filt)

        # Extract the pierce points
        trace_pierce = seis[0].stats['piercepoints'][
            'P'+str(depth)+'s'][str(depth)]
        
        #Make Shapely points
        point = Point(float(trace_pierce[1]), float(trace_pierce[2]))

        # Check whether event pierce points located within bounds of regions
        if box.contains(point):
            count_yes = count_yes + 1

            # Extract various bits of info
            epi_dist = seis[0].stats['dist']
            RF_amp = getattr(seis[0], filt)['iterativedeconvolution']

            # Find out which bin this should go into in the stack
            rounded_epi_dist = (np.floor(epi_dist))

            # Find the index in the stack matrix that this epicentral distance
            # relates to
            rounded_epi_dist_loc = int(np.floor((epi_dist-30)/step))

            # Add every RF amplitude (vector) to the correct column (bin) in
            # the STACK matrix
            if len(RF_amp) == 1751:
                try:
                    STACK[:, rounded_epi_dist_loc] = STACK[:, rounded_epi_dist_loc]+RF_amp
                except:
                    print('something wrong')
                    countwrong = countwrong+1
            else:
                print('not the right length!')

            # Add 1 to the count of this bin to keep track of how many RFs it
            # contains
            try:
                counter[rounded_epi_dist_loc] = counter[rounded_epi_dist_loc] + 1
            except:
                print('here as well')

        else:
            count_no = count_no + 1

# Print counters
print (count_yes)
print (count_no)


#----------------------------------------Plot the stack-------------------

# Get the time y axis (same in each RF)
time = getattr(seis[0], filt)['time']

# Used to emphasise the smaller peaks (other than direct P)
NORMALIZATION = 0.075

# Plot STACK
plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
# Set up first subplot
ax1 = plt.subplot(gs[1])

# Smoothing process
if smoothing:

    # Create new stack and counter
    STACK_NEW = np.zeros([1751, epi_steps-1])
    counter_new = np.zeros([epi_steps-1])

    # Loop through the stack and counter, careful of the ends of the arrays
    for p in range(len(counter)):
        if p == 0:
            STACK_NEW[:, p] = STACK[:, p]+STACK[:, p+1]
            counter_new[p] = counter[p]+counter[p+1]
        elif p == (len(counter) - 1):
            STACK_NEW[:, p] = STACK[:, p]+STACK[:, p-1]
            counter_new[p] = counter[p]+counter[p-1]
        else:
            STACK_NEW[:, p] = STACK[:, p]+STACK[:, p-1]+STACK[:, p+1]
            counter_new[p] = counter[p]+counter[p-1]+counter[p+1]

    # Normalization and average
    # Divide stack by number of RF in each bin - recorded in counter
    for m in range(len(counter)):
        if counter[m] > 0:
            STACK_NEW[:, m] = STACK_NEW[:, m]/counter_new[m]

            # Check that direct P has amp of 1
            STACK_NEW[:, m] = STACK_NEW[:, m]/(np.nanmax(np.abs(STACK_NEW[:, m])))

    # Normalize after direct P to bring out other features
    STACK_NEW[:,:] = STACK_NEW[:,:]/NORMALIZATION  

    # Plotting command (x axis, y axis, matrix of values, colourmap, min and
    # max values)
    plt.pcolor(epi_range, time, STACK_NEW, cmap='seismic', vmin=-1, vmax=1)

# No smoothing process
else:
    for m in range(len(counter)):
        if counter[m] > 0:
            STACK[:, m] = STACK[:, m]/counter[m]
            STACK[:, m] = STACK[:, m]/(np.nanmax(np.abs(STACK[:, m])))     
    STACK[:,:] = STACK[:,:]/NORMALIZATION  
    plt.pcolor(epi_range, time, STACK, cmap='seismic', vmin=-1, vmax=1)
    
# Axes labels
plt.ylabel('Time (seconds)', fontsize=24)
plt.xlabel('Epicentral Distance (degrees)', fontsize=24)

# Axes limits
plt.gca().set_xlim([30, 90])
plt.gca().set_ylim([0, 150])
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(np.arange(0,150,30), fontsize=20)


#----------------------------------------Plot Predicted Travel Times------

# Check if any curves to plot
if plot_travel_time_curves:

    # Split the input into individual phases
    phases = phases_to_plot

    # Loop through each of the phases to plot
    for t in range(len(phases)):
        print (phases[t])

        # Read in the predicted travel time differences from the appropriate
        # file
        data = np.genfromtxt(
    '../Tools/Moveout_with_epicentral_distance/TT_'+str(phases[t])+'.dat')

        # Get the epicentral distance, set as x, and the travel time
        # differences, set as y
        x = data[:, 0]
        y = data[:, 1]

        # Plot the curves on the stack, annotating each with the name of the
        # phase
        plt.plot(x, y, linewidth=2)
        plt.annotate(
    str(phases[t]),
     xy=(x[600],
     y[600]),
     xytext=(30,
     -5),
     textcoords='offset points')

# No curves to plot
else:
    pass


#----------------------------------------Plot epicentral distance histogra

# Plot histogram
if plot_histogram:

    # Add second subplot
    ax2 = plt.subplot(gs[0])

    # Loop through the bins
    for n in range(len(counter)):
        mindist = (n*step)+30
        maxdist = mindist+step

        # Plot rectangle of histogram ((startx,,starty),width, height))
        ax2.add_patch(patches.Rectangle((mindist, 0), step, counter[n]))

    # Axes labels
    plt.ylabel('Number of RFs', fontsize=24)
    plt.xlabel('Epicentral Distance (degrees)', fontsize=24)

    # Axes limits
    plt.gca().set_xlim([30, 90])
    plt.gca().set_ylim([0, np.nanmax(np.abs(counter[:]))])

    # Subplot arrangement
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    plt.subplots_adjust(hspace=0)
    
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

# Don't plot histogram
else:
    pass


#----------------------------------------Final plotting-------------------

# Set up a few aspects of the title
if smoothing:
    Smooth = 'S'
else:
    Smooth = ' '

# Title
plt.suptitle(
    'Epicentral Distance Stack - '+str(
        depth)+'km - No. of RFs: '+str(
            count_yes)+' \n Lat/Lon '+str(
                lat1)+'/'+str(
                    lat3)+'/'+str(
                        lon1)+'/'+str(
                            lon2)+' - Bin Size: '+str(
                                bin_size)+' - '+str(
                                    Smooth))

plt.savefig(savepath+'Epicentral_Distance_Stack_'+str(count_yes)+'_RFs'+'.png')

# Plot
plt.show()
