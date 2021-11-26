#

from matplotlib import lines
from obspy import read
import matplotlib.pyplot as plt
import time
import glob
import numpy as np
from obspy import UTCDateTime
import sys
import os.path

def manualremoval_perstation(Data, Filt):
# Command line help
    # if len(sys.argv) != 2 or str(sys.argv[1]).lower() == 'help':
    #     print('\n')
    #     print('-----------------------------------------------------------------------------------------------------------------------')
    #     print(sys.argv[0])
    #     print('-----------------------------------------------------------------------------------------------------------------------')
    #     print('Description:           [OPTIONAL] Manual removal of RFs per station from the selected_RFs.dat list via epicentral distance.')
    #     print('Inputs:                Data directory (usually ../Data/), filter band')
    #     print('Outputs:               n/a\n')
    #     print('Usage:                 >> python3 8_plot_data_perstation.py datadirectory filterband')
    #     print('Options [1]:           jgf1, jgf2, jgf3, tff1, tff2, tff3, tff4 or tff5 [str]')
    #     print('-----------------------------------------------------------------------------------------------------------------------')
    #     print('\n')
    #     sys.exit()

    direc = Data
    flag = 'SV'
    filt = str(Filt)
    counter = 0
    goodrflist = []
    stadirs = sorted(glob.glob(direc+'/*'))

    for stadir in stadirs:
        print(stadir, Filt, Model)

# ------------------------------------- ED --------------------------------------------- #

        valid = True
        while valid:

            print('EXAMINE PER STATION FOR RF VS EPICENTRAL DISTANCE')

            if os.path.getsize(stadir+'/selected_RFs.dat'):
                EDs = input("Input E.D. of RFs to remove:")

                sta = stadir.replace(direc+'/','')
                
                if os.path.isfile(stadir+'/selected_RFs_'+noisefilter+Filt+'.dat'):
                    if EDs == "all":
                        open(stadir+'/selected_RFs_'+noisefilter+Filt+'.dat', "w").close()
                        open(stadir+'/SNR_selected_RFs_'+noisefilter+Filt+'.dat', "w").close()
                        print('.dat files cleared.')
                        valid = False
                    elif EDs == "none":
                        print('.dat files unedited.')
                        valid = False
                    else:   
                        if ' ' in EDs:
                            EDs = EDs.split(' ')
                        elif ', ' in EDs:
                            EDs = EDs.split(', ')
                        elif ',' in EDs:
                            EDs = EDs.split(',')
                        else:
                            EDs = np.array(EDs)

                        print(EDs)
                
                        with open(stadir+'/selected_RFs_'+noisefilter+Filt+'.dat','r') as f:
                            for line in f:
                                for ED in EDs:
                                    if "_"+ED+"_" in line:
                                        line.replace(line + '\n', '')
                                        print('RF at '+ED+' degrees removed from .dat files.')
                        with open(stadir+'/SNR_selected_RFs_'+noisefilter+Filt+'.dat','r') as f:
                            for line in f:
                                for ED in EDs:
                                    if "_"+ED+"_" in line:
                                        line.replace(line + '\n', '')
                                        print('RF at '+ED+' degrees removed from SNR .dat files.')
                                        valid = False

                        goagain = input("Remove more? [Y/N]")
                        if goagain =='Y':
                            valid = True
                        elif goagain =='N':
                            valid = False

            else:
                print('No RFs to inspect.')
                valid = False
               
manualremoval_perstation(sys.argv[1], sys.argv[2])
