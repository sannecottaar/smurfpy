# AUTOMATED QUALITY CONTROL SCRIPT FOR RF ANALYSIS#################

# This script looks at computed RF and removes low quality ones based on
# set criteria:

# -if maximum value of the Receiver function is actually between 23 and 27 seconds from the P wave arrival (ie P wave is largest arrival)
# -if 60% of SV is recovered (when RF in reconvolved with Z)
# -if there are no strong peaks before and no huge peaks after the main P arrival
# -also computes SNR on *HZ, *HR - either ratio of P energy to trace (Sanne) or Pprior to Ppost energy (Alistair)
# -new SNR accepts waveforms where first P arrival is not highest amplitude.
# poor RF are moved to a sub folder and a list of automatically slected RF is produced:selected_RFs_jgf1.dat
#

# import modules

from obspy import read
from obspy.core import trace
import os.path
import glob
import numpy as np
import sys

# Command line help
if len(sys.argv) != 2 or str(sys.argv[1]).lower() == 'help':
    print('\n')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print(sys.argv[0])
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('Description:           Removes low quality ones based on set criteria:')
    print('                           1. Minimum percentage of radial compoment to be fit (after reconvolving the RF with the')
    print('                           vertical component (fitmin)')
    print('                           2. Peak amplitude max threshold before main P-wave arrival (noisebefore)')
    print('                           3. Peak amplitude max threshold after main P-wave arrival (noiseafter)')
    print('                           4. Peak amplitude min threshold after main P-wave arrival (minamp)')
    print('Inputs:                Data directory (usually ../Data/), horizontal component (usually radial), filter band, minimum/')
    print('                       maximum epicentral distances, SNR calculation type, fitmin, noisebefore, noiseafter, minamp')
    print('Outputs:               Two ".dat" files specific to the chosen filter band recording the good RF files and the good')
    print('                       RF file SNR ratios (V & R components)\n')
    print('Usage:                 >> python3 6_auto_select_receiver_functions.py filterband mindistance maxdistance')
    print('Options [1]:           jgf1, jgf2, jgf3, tff1, tff2, tff3, tff4 or tff5 [str]')
    print('Options [2,3]:         epicentral distance [int])
    print('Recommended:           python3 6_auto_select_receiver_functions.py jgf1 30 90')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('\n')
    sys.exit()

# negative length of file name in characters +1
#****
file_char_len = -42
#****
direc = '../Data'
flag = 'SV'
filt = str(sys.argv[1])
mindist = int(sys.argv[2])
maxdist = int(sys.argv[3])
count = 0

# set noise criteria
# Use Sanne Vert/Rad trace SNR measure
Sanne = False
# Use Alistair Vert/Rad trace SNR measure
Alistair = True
# Use Gao and Liu 2014 SNR calculation
Gao = False

    
fitmin = 60  # Minimum 60 percent of radial compoment should be fit (after reconvolving the RF with the vertical component
noisebefore = 0.4  # There should be no peaks before the main P wave which are more than 40% of the amplitude of the Pwave
noiseafter = 0.7  # There should be no peaks after the main P wave which are more than 70% of the amplitude of the Pwave
minamp = 0.04  # There should be signals after the P wave which are at least 4% of the main P wave. Otherwise something has gone wrong...

# Set counts
Noise = 0
LowE = 0
BeforeP = 0
AfterP = 0
Fit = 0
P_time = 0
Epi_dist = 0

stadirs = sorted(glob.glob(direc + '/*'))

for stadir in stadirs:
    # This is where selected data will be listed.
    goodrffile = open(stadir + '/selected_RFs_' + filt + '.dat', 'w')
    SNRrffile = open(stadir + '/SNR_selected_RFs_' + filt + '.dat', 'w')

    stalist = sorted(glob.glob(stadir + '/*.PICKLE'))

    # Make directories for rubbish data if it doesn't already exist
    rm_direc = stadir + '/auto_removed'
    if not os.path.exists(rm_direc):
        os.makedirs(rm_direc)

    c = 0
    # Loop through receiver functions
    for i in range(len(stalist)):
        print('checking quality', i + 1, "/", len(stalist), stalist[i])
        if os.path.isfile(stalist[i]):
            seis = read(stalist[i], format='PICKLE')
            dist = seis[0].stats['dist']
            if dist > mindist and dist < maxdist:  # Check distance range
                RF = trace.Trace()
  
                findRF = getattr(seis[0], filt)
                RF.data = np.real(findRF['iterativedeconvolution'])
                RF.data = RF.data / np.max(np.abs(RF.data))
                fitid = findRF['iterativedeconvolution_fit']
                # Test if maximum value of the Receiver function is actually
                # between 23 and 27 seconds from the P wave arriva
                indm = np.argmax(np.abs(RF.data))
                withinrange = True if (indm > 230 and indm < 270) else False
                if withinrange:
                    # If 60% of SV is recovered, and RF passes the peak test,
                    # save the RF
                    if fitid > fitmin:
                        # test if there are no strong peaks before and no huge
                        # peaks after the main arrival
                        if np.max(np.abs(RF.data[0:indm - 40])) < noisebefore:
                            if np.max(np.abs(RF.data[indm + 40:indm + 1200])) < noiseafter:
                                if np.max(np.abs(RF.data[indm + 200:-1])) > minamp:
                                    Ptime = seis[0].stats.traveltimes['P']  # set P arrival time
                                    Ptime = seis[0].stats.event.origins[0].time + Ptime
                                    fmax=findRF['maxfreq']
                                    fmin=findRF['minfreq']
                                    
                                    # Sanne SNR measure                                   
                                    if Sanne:
                                        vertical = seis.select(channel='*HZ')[0]
                                        Pref = vertical.slice(Ptime - 25., Ptime + 150.)  # Cut out P arrival on vertical
                                        radial = seis.select(channel='*HR')[0]
                                        SVref = radial.slice(Ptime - 25., Ptime + 150.)  # Cut out P arrival on radial
                                        Pref.filter('bandpass',freqmin=fmin,freqmax=fmax,corners=2,zerophase=True)
                                        SVref.filter('bandpass',freqmin=fmin,freqmax=fmax,corners=2,zerophase=True)
                                        Pref.taper(max_percentage=0.05, type='cosine')
                                        SVref.taper(max_percentage=0.05, type='cosine')
                                        
                                        # test mean energy around P arrival, vs.
                                        # rest of the trace
                                        Ponly = Pref.slice(Ptime - 5., Ptime + 20.)
                                        Penergy = np.sum(Ponly.data * Ponly.data) / float(len(Ponly.data))
                                        # print ('Penergy= ', Penergy)
                                        # Pref=Pref.slice(Ptime+20.,Ptime+90.)
                                        Tracenergy = np.sum(Pref.data * Pref.data) / float(len(Pref.data))
                                        # print ('Tracenergy= ', Tracenergy)
                                        Noisemeasure = Penergy / Tracenergy
                                        # test mean energy on radial component
                                        Tracenergy = np.sum(SVref.data * SVref.data) / float(len(SVref.data))
                                        Noisemeasure2 = Penergy / Tracenergy
                                        # set requirements for signal-to-noise
                                        # ratio
                                        if Noisemeasure > 2.5 and Noisemeasure2 > 2.0 and Noisemeasure2 < 100.:
                                            print(i,stalist[i], 'passed selection')
                                            count = count + 1
    
                                            print('good station')
                                            goodrffile.write(stalist[i] + '\n')  # List all RF filenames is good_RFs.dat
                                            SNRrffile.write(stalist[i] + Noisemeasure + Noisemeasure2 + '\n') # Print calculated SNR
                                            # if file doesn't meet any of these
                                            # criteria - move to thrown out folder
                                        else:
                                            print("Station: " + stadir + ", Event: " + stalist[i] + ", Poor signal to noise" + '\n')
                                            print("Noisemeasure = " + str(Noisemeasure) + ", Noisemeasure2 = " + str(Noisemeasure2) + '\n')
                                            Noise += 1
                                        
                                    # Alistair Noise measure 
                                    # RMS(Pre-arrival window)/RMS(post-arrival window)
                                    if Alistair:
                                        vertical = seis.select(channel='*HZ')[0]
                                        Pref = vertical.slice(Ptime - 100., Ptime + 150.)  # Cut out P arrival on vertical
                                        radial = seis.select(channel='*HR')[0]
                                        SVref = radial.slice(Ptime - 100., Ptime + 150.)  # Cut out P arrival on radial
                                        Pref.filter('bandpass',freqmin=fmin,freqmax=fmax,corners=2,zerophase=True)
                                        SVref.filter('bandpass',freqmin=fmin,freqmax=fmax,corners=2,zerophase=True)
                                        Pref.taper(max_percentage=0.05, type='cosine')
                                        SVref.taper(max_percentage=0.05, type='cosine')
                                        
                                        # Test mean energy prior to P v.s. After P. on vertical
                                        Pprior = Pref.slice(Ptime - 65., Ptime - 5.) # Select 1 min of data with safety gap
                                        Pprior_RMS = np.sqrt(np.sum(Pprior.data * Pprior.data) / float(len(Pprior.data)))                                        
                                        
                                        Ppost = Pref.slice(Ptime, Ptime + 60.) # Select 1 min of data with safety gap
                                        Ppost_RMS = np.sqrt(np.sum(Ppost.data * Ppost.data) / float(len(Ppost.data)))                                      
                                        
                                        SNR_V = Ppost_RMS / Pprior_RMS
                                        
                                        # Test mean energy prior to P v.s. After P. on vertical
                                        Pprior = SVref.slice(Ptime - 65., Ptime - 5.) # Select 1 min of data with safety gap
                                        Pprior_RMS = np.sqrt(np.sum(Pprior.data * Pprior.data) / float(len(Pprior.data)))                                        
                                        
                                        Ppost = SVref.slice(Ptime, Ptime + 60.) # Select 1 min of data with safety gap
                                        Ppost_RMS = np.sqrt(np.sum(Ppost.data * Ppost.data) / float(len(Ppost.data)))                                      
                                        
                                        SNR_R = Ppost_RMS / Pprior_RMS                                        
                                        
                                        # set requirements for signal-to-noise
                                        # ratio THink 2.5/1.75 works well for African data
                                        if SNR_V >= 2.5 and SNR_R >= 1.75:
                                            print(i,stalist[i], 'passed selection')
                                            count = count + 1
    
                                            print('good station')
                                            goodrffile.write(stalist[i] + '\n')  # List all RF filenames is good_RFs.dat
                                            SNRrffile.write(stalist[i] + ' ' + str(round(SNR_V,2)) + ' ' + str(round(SNR_R,2)) + '\n') # Print calculated SNR
                                            # if file doesn't meet any of these
                                            # criteria - move to thrown out folder
                                        else:
                                            print("Station: " + stadir + ", Event: " + stalist[i] + ", Poor signal to noise" + '\n')
                                            print("SNR_V = " + str(SNR_V) + ", SNR_R = " + str(SNR_R) + '\n')
                                            Noise += 1
                                            
                                            
                                            
                                    # Gao and Liu 2014 (JGR) Noise measure
                                    # SNR_V > 4
                                    # max(Pdata[-8s:17s])/mean(abs(Pdata[-20s:-10s]))
                                    if Gao:
                                        vertical = seis.select(channel='*HZ')[0]
                                        Pref = vertical.slice(Ptime - 100., Ptime + 150.)  # Cut out P arrival on vertical
                                        Pref.filter('bandpass',freqmin=fmin,freqmax=fmax,corners=2,zerophase=True)
                                        Pref.taper(max_percentage=0.05, type='cosine')
                                        
                                        #Find max of Vertical trace around Ptime
                                        Ppeak = Pref.slice(Ptime - 8., Ptime + 17.) # Select max amp predicted time
                                        Ppeak_max = np.max(Ppeak.data)
                                        # Get mean absolute amp prior to arrival
                                        Pprior = Pref.slice(Ptime - 20., Ptime - 10.)
                                        Pprior_mean = np.sum(np.abs(Pprior.data)) / float(len(Pprior.data)) 
                                        
                                        SNR_V = Ppeak_max / Pprior_mean
                                        SNR_R = 1
                                        # set requirements for signal-to-noise
                                        # ratio THink 2.5/1.85 works well.
                                        if SNR_V >= 4 and SNR_R >= 1:
                                            print(i,stalist[i], 'passed selection')
                                            count = count + 1
    
                                            print('good station')
                                            goodrffile.write(stalist[i] + '\n')  # List all RF filenames is good_RFs.dat
                                            SNRrffile.write(stalist[i] + ' ' + str(round(SNR_V,2)) + ' ' + str(round(SNR_R,2)) + '\n') # Print calculated SNR
                                            # if file doesn't meet any of these
                                            # criteria - move to thrown out folder
                                        else:
                                            print("Station: " + stadir + ", Event: " + stalist[i] + ", Poor signal to noise" + '\n')
                                            print("SNR_V = " + str(SNR_V) + ", SNR_R = " + str(SNR_R) + '\n')
                                            Noise += 1 
                                else:
                                    print("Station: " + stadir + ", Event: " + stalist[i] + ", Peaks < 4% after P-wave" + '\n')
                                    print("Maximum amplitude after P: " + str(np.max(np.abs(RF.data[indm + 200:-1]))) + '\n')
                                    LowE += 1
                            else:
                                print("Station: " + stadir + ", Event: " + stalist[i] + ", Peaks > 70% after P-wave" + '\n')
                                print("Maximum amplitude after P: " + str(np.max(np.abs(RF.data[indm + 40:indm + 1200]))) + '\n')
                                AfterP += 1
                        else:
                            print("Station: " + stadir + ", Event: " + stalist[i] + ", Peaks > 40% before P-wave" + '\n')
                            print("Maximum amplitude before P: " + str(np.max(np.abs(RF.data[0:indm - 40]))) + '\n')
                            BeforeP +=1
                    else:
                        print("Station: " + stadir + ", Event: " + stalist[i] + ", Less than 60% fit" + '\n')
                        print("Percentage fit: " + str(fitid) + '\n')
                        Fit += 1                        
                else:
                    print("Station: " + stadir + ", Event: " + stalist[i] + ", P-wave not between 23 and 27 seconds" + '\n')
                    print("P arrival-time: " + str(indm) + '\n')
                    P_time +=1
            else:
                print("Station: " + stadir + ", Event: " + stalist[i] + ", Not between " + str(mindist) + "-" + str(maxdist) + "deg " + '\n')
                print("Dist: " + str(dist) + '\n')
                Epi_dist +=1                    
    goodrffile.close()
    SNRrffile.close()
# Write counts and close
Total_count = Noise + LowE + AfterP + BeforeP + Fit + P_time + Epi_dist
print("\\\\\    Total new good RFs is " + str(count) + '\n')
print("Total RFs rejected: " + str(Total_count) + '\n')
print("Noise: " + str(Noise) + " Low Energy: " + str(LowE) + " Amplitude after P: " + str(AfterP) + '\n')
print("Amplitude before P: " + str(BeforeP) + " Percentage fit: " + str(Fit) + " P arrival-time: " + str(P_time) + '\n')
print("Epicentral Distance: " + str(Epi_dist) + '\n')
