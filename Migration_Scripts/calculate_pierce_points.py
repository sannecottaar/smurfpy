'''
Piercepoints_Calc.py
This script computes predicted pierce points for user defined phases based on TauP.
It loops through all the data, checking whether the pierce points have been calculated,
and if not working them out, writing the data to sperate file and the PICKLE file itself.

Phases can be computed for every discontinuity depth and evaluated at every discontinuity
depth in prem_added_discon_taup.npz or ak135_added_discon_taup.npz.

sys.argv[1] = Depth of piercepoints
sys.argv[2] = Phase
sys.argv[3] = Filter

eg. Piercepoints_Calc.py 410 P410s jgf1
'''

# Import all the relevant modules
from obspy import read
import os, sys, glob
from obspy.taup import TauPyModel

# Command line help
if len(sys.argv) != 4 or str(sys.argv[1]).lower() == 'help':
    print('\n')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print(sys.argv[0])
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('Description:           Calculate converted phase pierce points at discontinuity depths')
    print('Inputs:                Depth of piercepoints, Phase, filter band, 1D velocity model')
    print('Outputs:               Adds PP for given phase and discont depth to each Pickle file, prints to file')
    print('                       PP_DEPTHkm_PHASE_FILTER.txt\n')
    print('Usage:                 python3 calculate_pierce_points.py depth phase filter')
    print('Format [1]:            depth (km) [int]')
    print('Format [2]:            seismic phase [str]')
    print('Options [3]:           jgf1, jgf2, jgf3, tff1, tff2, tff3, tff4 or tff5 [str]')
    print('Recommended:           python3 calculate_pierce_points.py 410 P410s jgf1')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('\n')
    sys.exit()
    

PREM=True
ak135=False
if PREM:
    mod='prem_added_discon_taup.npz'
elif ak135:
    mod='ak135_added_discon_taup.npz'
    
# Allows for more depths than normal PREM
# Discrepancy only appears in 4th decimal of lat/lon for P400s
taupmodel = TauPyModel(model=mod)

# Make and open a data file in the Piercepoints directory for this phase
# and depth ('w' = writable)
PP = open('PP_' + str(sys.argv[1]) + 'km_' + str(sys.argv[2]) + '_' + str(sys.argv[3]) + '_test.txt', 'w')

# Loop through stations
direc = '../RF_Data'
stadirs = glob.glob(direc + '/IU*')
for stadir in stadirs:
    stalist = []
    if os.path.isfile(stadir + '/selected_RFs_'+str(sys.argv[3])+'.dat'):
        with open(stadir + '/selected_RFs_'+str(sys.argv[3])+'.dat') as a:
            starfs = a.read().splitlines()
            for line in starfs:
                stalist.append(line)

    # Loop through events
    for s in range(len(stalist)):
        seis = read(stalist[s], format='PICKLE')
        print (stalist[s])

        # Note stats of stream
        EVLA = seis[0].stats['evla']
        EVLO = seis[0].stats['evlo']
        EVDP = seis[0].stats['evdp']
        STLA = seis[0].stats['stla']
        STLO = seis[0].stats['stlo']

        # Check if there's a piercepoints dictionary
        if not hasattr(seis[0].stats, 'piercepoints'):

            # Start a piercepoints dictionary
            seis[0].stats.piercepoints = dict()

        # Set up a variable to check whether the process needs to be done
        runthis = False

        # Check if the dictionary has this phase in
        if not str(sys.argv[2]) in seis[0].stats.piercepoints.keys():
            runthis = True

        # Check whether the dictionary with this phase has this depth in
        if str(sys.argv[2]) in seis[0].stats.piercepoints.keys():
            if not str(sys.argv[1]) in seis[0].stats.piercepoints[sys.argv[2]].keys():
                runthis = True

        # If the file neither has the phase or depth, run the rest of the
        # script
        if runthis:
            print('runthis is True')
            # Call taup_pierce to get piercepoints
            arrivals = taupmodel.get_pierce_points_geo(EVDP,EVLA,EVLO,STLA,STLO,phase_list=([sys.argv[2]]))

            # Index from source side to find nearest pierce depth
            pierce_list = arrivals[0].pierce[::-1]
            # Find nearest depth
            index = min(range(len(pierce_list)), key=lambda i: abs(float(pierce_list[i][3])-float(sys.argv[1])))

            # Make a separate dictionary within piercepoints for each phase
            seis[0].stats.piercepoints[sys.argv[2]] = dict()

            # For each phase, create a dictionary for each depth, writing
            # piercepoint depth/lat/lon info
            seis[0].stats.piercepoints[sys.argv[2]][
                sys.argv[1]] = [pierce_list[index][3],pierce_list[index][4],pierce_list[index][5]]

            # Write pierce point data to original pickle file
            seis.write(stalist[s], format='PICKLE')

            # Add line to opened output file followed by new line
            PP.write(str(float(pierce_list[index][3])) + ' ' + str(round(pierce_list[index][4],2)) + ' '
                     + str(round(pierce_list[index][5],2)) + '\n')

        else:
            PP.write(seis[0].stats.piercepoints[sys.argv[2]][sys.argv[1]][0] + ' ' + \
            seis[0].stats.piercepoints[sys.argv[2]][sys.argv[1]][1] + ' ' + \
            seis[0].stats.piercepoints[sys.argv[2]][sys.argv[1]][2] + '\n')

# Close the file you are writing to
PP.close()
print ('Finished')
