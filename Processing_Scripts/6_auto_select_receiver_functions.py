##### AUTOMATED QUALITY CONTROL SCRIPT FOR RF ANALYSIS#################

# This script looks at computed RF and removes low quality ones based on set criteria:
 
# -if maximum value of the Receiver function is actually between 24 and 26 seconds from the P wave arrival (ie P wave is largest arrival)
# -if 60% of SV is recovered (when RF in reconvolved with Z)
# -if there are no strong peaks before and no huge peaks after the main P arrival

# poor RF are moved to a sub folder and a list of automatically slected RF is produced:selected_RFs_jgf1.dat
##################################################################

# import modules

from obspy import read
from obspy.core import trace
import os.path
import glob
import shutil
import numpy as np


#negative length of file name in characters +1
#****
file_char_len=-51
#****
direc = 'DataRF/'
flag = 'SV'
filt= 'jgf1'
count=0

fitmin = 60 # Minimum 60 percent of radial compoment should be fit (after reconvolving the RF with the vertical component
noisebefore = 0.3 # There should be no peaks before the main P wave which are more than 30% of the amplitude of the Pwave
noiseafter = 0.7# There should be no peaks after the main P wave which are more than 70% of the amplitude of the Pwave
minamp = 0.04 # There should be signals after the P wave which are at least 4% of the main P wave. Otherwise something has gone wrong...


stadirs= glob.glob(direc+'/*')
for stadir in stadirs:
    # This is where selected data will be listed. 
    goodrffile= open(stadir+'/selected_RFs_jgf1.dat','w')

    stalist=glob.glob(stadir+'/*.PICKLE') 

    #Make directories for rubbish data if it doesn't already exist
    rm_direc= stadir+'/auto_removed'
    if not os.path.exists(rm_direc):
        os.makedirs(rm_direc)

    c=0
    # Loop through receiver functions
    for i in range(len(stalist)):
        print('checking quality', i+1,"/",len(stalist), stalist[i])
        if os.path.isfile(stalist[i]):
            seis=read(stalist[i],format='PICKLE')
            dist=seis[0].stats['dist']
            if dist>30 and dist <90: # Check distance range
                RF=trace.Trace()
                #try:
                findRF = getattr(seis[0],filt)
                RF.data= np.real(findRF['iterativedeconvolution'])
                RF.data=RF.data/np.max(np.abs(RF.data))
                fitid=findRF['iterativedeconvolution_fit']
                print(fitid)
                # Test if maximum value of the Receiver function is actually between 24 and 26 seconds from the P wave arrival
                indm=np.argmax(np.abs(RF.data))
                print(indm)
                withinrange= True if (indm>230 and indm <270) else False
                if withinrange:
                    # If 60% of SV is recovered, and RF passes the peak test, save the RF
                    if fitid>fitmin:
                        # test if there are no strong peaks before and no huge peaks after the main arrival
                        if np.max(np.abs(RF.data[0:indm-40]))<noisebefore:
                         if np.max(np.abs(RF.data[indm+40:indm+1200]))<noiseafter:
                          if np.max(np.abs(RF.data[indm+200:-1]))>minamp:
                              Ptime=seis[0].stats.traveltimes['P'] # set P arrival time
                              Ptime=seis[0].stats['starttime']+Ptime
                              vertical = seis.select(channel='BHZ')[0]
                              Pref=vertical.slice(Ptime-25.,Ptime+150.) # Cut out P arrival on vertical
                              radial = seis.select(channel='*HR')[0]
                              SVref=radial.slice(Ptime-25.,Ptime+150.) # Cut out P arrival on radial

                              # test mean energy around P arrival, vs. rest of the trace
                              Ponly=Pref.slice(Ptime-5.,Ptime+20.)
                              Penergy=np.sum(Ponly.data*Ponly.data)/float(len(Ponly.data))
                              #print ('Penergy= ', Penergy)                              
                              #Pref=Pref.slice(Ptime+20.,Ptime+90.)
                              Tracenergy=np.sum(Pref.data*Pref.data)/float(len(Pref.data))
                              #print ('Tracenergy= ', Tracenergy)                        
                              Noisemeasure=Penergy/Tracenergy
                              # test mean energy on radial component
                              Tracenergy=np.sum(SVref.data*SVref.data)/float(len(SVref.data))
                              Noisemeasure2=Penergy/Tracenergy
                              # set requirements for signal-to-noise ratio
                              if Noisemeasure > 2.5 and Noisemeasure2 >2.0 and Noisemeasure2 <100. :
                                  print stalist[i], 'passed selection'
                                  count=count+1
                              
                                  print('good station')
                                  goodrffile.write(stalist[i]+'\n') # List all RF filenames is good_RFs.dat
                                  #if file doesn't meet any of these criteria - move to thrown out folder
                              else:
                              shutil.move(stalist[i],rm_direc+stalist[i][file_char_len:])
                          else:
                            shutil.move(stalist[i],rm_direc+stalist[i][file_char_len:])
                         else:
                          shutil.move(stalist[i],rm_direc+stalist[i][file_char_len:])
                        else:
                          shutil.move(stalist[i],rm_direc+stalist[i][file_char_len:])
                    else:
                      shutil.move(stalist[i],rm_direc+stalist[i][file_char_len:])
                else:
                  shutil.move(stalist[i],rm_direc+stalist[i][file_char_len:])
  
                 #except:
                 #  print('test failed for' + stalist[i])

    goodrffile.close()


