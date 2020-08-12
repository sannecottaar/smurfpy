# RAW SEISMOGRAM PLOTTING SCRIPT FOR RF ANALYSIS#################

# This script plots ZRT seismogram components with predicted arrival times
# after initial processing and before RF creation


from obspy import read
import matplotlib.pyplot as plt
import glob
import numpy as np

direc = '../Data'
stadirs = glob.glob(direc + '/*')

for stadir in stadirs:

    # Loop through events
    stalist = glob.glob(stadir + '/*.PICKLE')

    # Loop through data
    if(len(stalist) > 0):
        for i in range(0, len(stalist),20):  # Plotting every 20
                print(i)
                seis = read(stalist[i], format='PICKLE')

                # Extract epi dist of event/station from header info
                distdg = seis[0].stats['dist']

                #----------------PLOT VERTICAL--------------------------------
                # 1 of 3 subplots (grid 1 row 3 columns)
                plt.subplot(1, 3, 1)
                vertical = seis.select(channel='*HZ')[
                    0]  # Extract V from stream
                vertical.filter(
                    'bandpass',
                    freqmin=0.01,
                    freqmax=.1,
                    corners=2,
                    zerophase=True)

                plt.plot(
                    vertical.times(),
                    vertical.data / np.max(vertical.data) + np.round(distdg),
                    'k')  # Plot normalised data as black line

                # Plot predicted travel times as coloured points (add more
                # phases as wanted)
                plt.plot(
                    seis[0].stats.traveltimes['P'],
                    np.round(distdg),
                    '.b')
                plt.plot(
                    seis[0].stats.traveltimes['S'],
                    np.round(distdg),
                    '.g')

                # Label plot
                plt.title('Vertical')
                plt.ylabel('Epi dist (degrees)')
                plt.xlabel('Time from Origin (sec)')

                #----------------PLOT RADIAL--------------------------------
                # 2 of 3 subplots
                plt.subplot(1, 3, 2)
                radial = seis.select(channel='*HR')[0]

                radial.filter(
                    'bandpass',
                    freqmin=0.01,
                    freqmax=.1,
                    corners=2,
                    zerophase=True)

                plt.plot(
                    radial.times(),
                    radial.data / np.max(radial.data) + np.round(distdg),
                    'k')

                plt.plot(
                    seis[0].stats.traveltimes['P'],
                    np.round(distdg),
                    '.b')
                plt.plot(
                    seis[0].stats.traveltimes['P410s'],
                    np.round(distdg),
                    '.b')
                plt.plot(
                    seis[0].stats.traveltimes['P660s'],
                    np.round(distdg),
                    '.b')
                plt.plot(
                    seis[0].stats.traveltimes['S'],
                    np.round(distdg),
                    '.g')

                plt.title('Radial')
                plt.xlabel('Time from Origin (sec)')

                #----------------PLOT TRANSVERSE--------------------------------
                # 3 of 3 subplots
                plt.subplot(1, 3, 3)
                transverse = seis.select(channel='*HT')[0]
                transverse.filter(
                    'bandpass',
                    freqmin=0.01,
                    freqmax=.1,
                    corners=2,
                    zerophase=True)

                plt.plot(
                    transverse.times(),
                    transverse.data /
                        np.max(transverse.data) + np.round(distdg),
                    'k')
                plt.plot(
                    seis[0].stats.traveltimes['P'],
                    np.round(distdg),
                    '.b')
                plt.plot(
                    seis[0].stats.traveltimes['S'],
                    np.round(distdg),
                    '.g')
                plt.title('Transverse')
                plt.xlabel('Time from Origin (sec)')

        plt.show()
