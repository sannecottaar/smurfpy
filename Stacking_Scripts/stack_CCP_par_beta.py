# common conversion point stacking

# Uses routines from common_converion_point_stack.
# Needs a directory '../CCP_Volumes'

#----------------------------------------------------#
import common_conversion_point_stack_par_beta as CCP_stack
import os, glob
import time
import sys
from obspy import read
import concurrent.futures
import numpy as np
import msgpack
import msgpack_numpy as m
m.patch()

#----------------------------------------------------#

# Command line help
if len(sys.argv) != 11 or str(sys.argv[1]).lower() == 'help':
    print('\n')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print(sys.argv[0])
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('Description:           Wrapper for the function contained within common_conversion_point_stack_par_beta.py')
    print('Inputs:                depth conversion, lon/lat box, filter band')
    print('Outputs:               Common conversion point stack volume (PICKLE)\n')
    print('Usage:                 >> python3 stack_CCP_par_beta.py conversion lonmin lonmax latmin latmax rffilter newstack cores')
    print('Format                 1,2: [str], 3-6: [float], 7: [str], 8: [float], 9: [bool], 10: [int]')
    print('Recommended:           >> python3 stack_CCP_par_beta.py CCP_Global prem -179.0 179.0 -89.0 89.0 jgf1 2.0 True 12')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('\n')
    sys.exit()

# Initial options
name = str(sys.argv[1])
conversion = str(sys.argv[2]) # conversion use
lonmin = float(sys.argv[3]) 
lonmax = float(sys.argv[4]) 
latmin = float(sys.argv[5]) 
latmax = float(sys.argv[6])
rffilter=str(sys.argv[7])  # RF filter used
factor = float(sys.argv[8])  # e.g. 2.0 to double the fresnel zone width
newstack = bool(str(sys.argv[9])) #Starting a new stack (True) or adding data to an old one (False)
cores = int(sys.argv[10])

depmin=60.
depmax=1300.
lonrez=(lonmax-lonmin)*2.+1. # every 0.5 degrees
latrez=(latmax-latmin)*2.+1. # every 0.5 degrees
deprez=(depmax-depmin)/2. # every 2 km

script_start = time.time()

## Intialize stack or load latest stack
CCP_master= CCP_stack.ccp_volume()#

CCP_master.start_empty_volume_master(name= name, filter=rffilter, conversion=conversion, factor=factor, lonmin=lonmin, lonmax=lonmax, lonrez=lonrez, latmin=latmin, latmax=latmax, latrez=latrez,depmin=depmin, depmax=depmax, deprez=deprez)

stations=glob.glob('../Data/*')

print('Starting loop to add RFs...')
for idx, sta in enumerate(stations):
    print(sta)
    rflist=[]
    if os.path.exists(sta + '/selected_RFs_'+str(rffilter)+'.dat'):
        station_file = open(sta + '/selected_RFs_'+str(rffilter)+'.dat','r')
        file_list=station_file.read().strip().split()
        for file in file_list:
            seis=read(file, format='PICKLE')
            if hasattr(seis[0],rffilter):
                rflist.append(file)
        
        def mk_sub_stack(index):
            global rflist, name, rffilter, conversion, factor, idx,lonmin, lonmax, lonrez, latmin, latmax, latrez, depmin, depmax, deprez
            CCP_sub=CCP_stack.ccp_volume()
            outfilename='../CCP_volumes/'+name+'_'+rffilter+'_'+conversion+'_'+str(factor)+'/Stack_'+str(int(idx))+'_'+str(int(index))+'.PICKLE'
            CCP_sub.start_empty_sub_volume(name=name, filter=rffilter, conversion=conversion, factor=factor, sta_num=int(idx), file_num=int(index), lonmin=lonmin, lonmax=lonmax, lonrez=lonrez, latmin=latmin, latmax=latmax, latrez=latrez,depmin=depmin, depmax=depmax, deprez=deprez)
            # Want to add the rflist[index] to the CCP_sub
            rffile=rflist[index]
            CCP_sub.addlist_sub(rffile=rffile,name=name,filter=rffilter,conversion=conversion,factor=factor, sta_num=int(idx), file_num=int(index))
            return index, rffile, outfilename

        def mk_master_stack(result):
            global CCP_master, name, rffilter, conversion, factor, idx, lonmin, lonmax, lonrez, latmin, latmax, latrez, depmin, depmax, deprez    
            print('adding: '+str(result))+' to CCP_master'
            # Record the added RF files in the RF_lists files.
            rffile=open('../CCP_volumes/'+name+'_'+rffilter+'_'+conversion+'_'+str(factor)+'/RF_lists/rflist'+(str(idx))+'.dat','a')
            rffile.write("%d %s \n" % (int(result[0]),result[1]))
            rffile.close()
            # Record the added Stack files in filenames.dat
            stackfile=open('../CCP_volumes/'+name+'_'+rffilter+'_'+conversion+'_'+str(factor)+'/filenames_subvol.dat','a')
            stackfile.write("%s \n" % (result[2]))
            stackfile.close()
            CCP_sub=CCP_stack.ccp_volume()
            # read in latest sub volume
            CCP_sub.load_latest_sub(name=name,filter=rffilter,conversion=conversion,factor=factor, sta_num=int(idx), file_num=int(result[0]))
            # read in latest master volume
            CCP_master.load_latest_master(name=name,filter=rffilter,conversion=conversion,factor=factor)

            # CCP_master.VOL['num']          = CCP_master.VOL['num']            + CCP_sub.VOL['num']             # count number of receiver functions
            CCP_master.VOL['volume']       = CCP_master.VOL['volume']         + CCP_sub.VOL['volume']          # stack receiver function into volume  
            CCP_master.VOL['volumeweight'] = CCP_master.VOL['volumeweight']   + CCP_sub.VOL['volumeweight']    # stack weights
            CCP_master.VOL['weightedvolumesquares'] = CCP_master.VOL['weightedvolumesquares']   + CCP_sub.VOL['weightedvolumesquares'] # Summed weighted volume squares
            CCP_master.VOL['volumesign']   = CCP_master.VOL['volumesign']     + CCP_sub.VOL['volumesign']      # stack sign of receiver function

            outfilename='../CCP_volumes/'+name+'_'+rffilter+'_'+conversion+'_'+str(factor)+'/Stack_master.PICKLE'
            with open(outfilename,'wb') as handle:
                msgpack.pack(CCP_master.VOL,handle,use_bin_type=True)
                handle.close()


        def del_sub_stack(result):
            global name, rffilter, conversion, factor, lonmin, lonmax, lonrez, latmin, latmax, latrez, depmin, depmax, deprez    
            print('Finished add sub volume... now cleaning...')
            try:
                if os.path.isfile(result[2]):
                    print('deleting '+str(result[2])+'...')
                    os.remove(result[2])
            except:
                print('Failed to delete '+str(result[2]))           


        index = range(0,len(rflist),1)
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            results = executor.map(mk_sub_stack, index)
            # Loop over each result and add to developing depth mask
            for result in results:
                mk_master_stack(result)
                del_sub_stack(result)

        def calc_master_errors(name='Megavolume',filter='jgf1',conversion='ak135',factor=2.):
            # Want to make the volumesigma for the CCP master here.
            # Defined as [Sum of weighted squared values] - [Sum of weighted values]**2 / [Sum of weights]
            CCP_master.load_latest_master(name=name,filter=rffilter,conversion=conversion,factor=factor)
            CCP_master.VOL['volumesigma']= CCP_master.VOL['weightedvolumesquares'] - ((CCP_master.VOL['volume']**2)/ CCP_master.VOL['volumeweight'])
            outfilename='../CCP_volumes/'+name+'_'+rffilter+'_'+conversion+'_'+str(factor)+'/Stack_master.PICKLE'
            with open(outfilename,'wb') as handle:
                msgpack.pack(CCP_master.VOL,handle,use_bin_type=True)
                handle.close()

        print('Calculating volumesigma for:',name,rffilter,conversion,factor)
        calc_master_errors(name=name,filter=rffilter,conversion=conversion,factor=factor)


print('Parallel CCP stacking for '+str(name)+'_'+str(rffilter)+'_'+str(conversion)+'_'+str(factor)+' is COMPLETE!!!')
print('It took', time.time()-script_start, 'seconds.')









