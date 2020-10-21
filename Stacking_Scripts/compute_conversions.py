import obspy
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

def compute_arrivals(conversion='prem'):

    if conversion=='prem':
        taupmodel = obspy.taup.taup_create.TauPCreate('../Tools/MODELS/PREM_FILES/prem.nd', output_filename='onediscon')
    elif conversion=='ak135':
        taupmodel = obspy.taup.taup_create.TauPCreate('../Tools/MODELS/ak135_FILES/ak135_added_discon_taup.tvel', output_filename='onediscon')
    
    vmodel= taupmodel.load_velocity_model()
    # introduce very minor disconinuities throughout the model 
    for layer in vmodel.layers:

        if layer[0]>0 and layer[0]<1500:
            layer['bot_p_velocity']=layer['bot_p_velocity']-0.00001
    discs=vmodel.get_discontinuity_depths()
    # Write out model
    taupmodel.run()

    # Start model in format to compute travel times from
    taupmodel =TauPyModel(model='./alldiscon.npz')
    # Compute reference phase
    ref =  taupmodel.get_travel_times(10,60,['P'])

    # computer travel times and slownesses
    directconversionT= [0]
    directconversionS =[0]    
    PmultipleT = [0]
    PmultipleS = [0]
    SmultipleT = [0]
    SmultipleS =[0]
    depthsdc = [0]
    depthsmp =[0]
    depthsms =[0]
    for d in discs:
        
        arr = taupmodel.get_travel_times(10,60,['P'+str(int(d))+'s'])
#        print (d,len(arr))
        if len(arr)>0:
            arr2 = taupmodel.get_travel_times(10,60,['P^0Pv'+str(int(d))+'p'])
            arr3 = taupmodel.get_travel_times(10,60,['P^0Pv'+str(int(d))+'s'])                      

            depthsdc.append(d)
            directconversionT.append(arr[0].time-ref[0].time)
            directconversionS.append((arr[0].ray_param-ref[0].ray_param)*np.pi/180.)           

            if len(arr2)>0:
                depthsmp.append(d)
                PmultipleT.append(arr2[0].time - ref[0].time)
                PmultipleS.append((arr2[0].ray_param-ref[0].ray_param)*np.pi/180.)
                depthsms.append(d)              
                SmultipleT.append(arr3[0].time - ref[0].time)
                SmultipleS.append((arr3[0].ray_param-ref[0].ray_param)*np.pi/180.)  
       
    depthsdc = np.array(depthsdc)
    depthsmp = np.array(depthsmp)
    depthsms = np.array(depthsms)
    directconversionT = np.array(directconversionT)
    directconversionS = np.array(directconversionS)
    PmultipleT  = np.array(PmultipleT)
    PmultipleS = np.array(PmultipleS)
    SmultipleT  = np.array(SmultipleT)
    SmultipleS = np.array(SmultipleS)
    


    # save files
    np.savetxt('time_slowness_directconversion.txt',list(zip(depthsdc,directconversionT, directconversionS)),fmt='%.6f',delimiter ='\t')
    np.savetxt('time_slowness_Pmultiple.txt',list(zip(depthsmp,PmultipleT, PmultipleS)),fmt='%.6f',delimiter ='\t')
    np.savetxt('time_slowness_Smultiple.txt',list(zip(depthsms,SmultipleT, SmultipleS)),fmt='%.6f',delimiter ='\t')   

    return((depthsdc,directconversionT, directconversionS), (depthsmp,PmultipleT, PmultipleS), (depthsms,SmultipleT, SmultipleS))


def compute_arrivals_one_depth(depth,conversion='prem'):

    if conversion=='prem':
        taupmodel = obspy.taup.taup_create.TauPCreate('../Tools/MODELS/PREM_FILES/prem.nd', output_filename='onediscon')
        vmodel= taupmodel.load_velocity_model()
        
    elif conversion=='ak135':
        taupmodel = obspy.taup.taup_create.TauPCreate('../Tools/MODELS/ak135_FILES/ak135_added_discon_taup.tvel', output_filename='onediscon')
        vmodel= taupmodel.load_velocity_model()



    # introduce one discontinuity at random depth
    # somewhat inspired by script from Github issue #1729 by Jack Walpole
    layers = np.array(vmodel.layers)

    rowdeeper = np.argmax(layers['top_depth']>=depth)

    newlayer =[depth]
    for i in [2,4,6,8,10]:
        f = scipy.interpolate.interp1d(layers['top_depth'],layers[layers.dtype.names[i]],kind='cubic')
        newlayer.append(f(depth))

    # Change top of layer to specific depth
    layers['top_depth'][rowdeeper]=depth
    c= 0
    for i in [2,4,6,8,10]:
            c=c+1
            layers[layers.dtype.names[i]][rowdeeper]=newlayer[c]
    print(layers[rowdeeper])

    # Change bottom of layer above to specific depth
    layers['bot_depth'][rowdeeper-1]=depth
    c= 0
    e = 0.000001 # shift to P wave velocity to create discontinuity
    for i in [3,5,7,9,11]:
            c=c+1
            layers[layers.dtype.names[i]][rowdeeper-1]=newlayer[c]
            if i ==3:
                layers[layers.dtype.names[i]][rowdeeper-1]=layers[layers.dtype.names[i]][rowdeeper-1]-e


    vmodel.layers = layers

    discs=vmodel.get_discontinuity_depths()

    # Build and reload model
    taupmodel.run()
    taupmodel =TauPyModel(model='./onediscon.npz')
    discs=taupmodel.model.s_mod.v_mod.get_discontinuity_depths()

    # Compute arrivals
    ref =  taupmodel.get_travel_times(10,60,['P'])
    arr = taupmodel.get_travel_times(10,60,['P'+str(int(depth))+'s'])
    arr2 = taupmodel.get_travel_times(10,60,['P^0Pv'+str(int(depth))+'p'])
    arr3 = taupmodel.get_travel_times(10,60,['P^0Pv'+str(int(depth))+'s'])

    # Return referenc models
    try:
        return( (depth, arr[0].time-ref[0].time, (arr[0].ray_param-ref[0].ray_param)*np.pi/180.),(depth, arr2[0].time-ref[0].time, (arr2[0].ray_param-ref[0].ray_param)*np.pi/180.),(depth, arr3[0].time-ref[0].time,(arr3[0].ray_param-ref[0].ray_param)*np.pi/180.))
    except:
        return( (depth, arr[0].time-ref[0].time, (arr[0].ray_param-ref[0].ray_param)*np.pi/180.))
        


if __name__=="__main__":

    # plot direct conversions and multiples
    
    lines= compute_arrivals()
    plt.plot(lines[0][1] ,lines[0][2],'k')   
    plt.plot(lines[1][1], lines[1][2],'r')
    plt.plot(lines[2][1], lines[2][2],'b')
    plt.xlim([-2,140])
    plt.ylim([-0.7,0.7])
    
    # plot specific depths
    for d in [340,445,650]:
        arr = compute_arrivals_one_depth(d)
        
        plt.plot(arr[0][1] , arr[0][2],'.g', markersize=10)
        plt.text(arr[0][1] , arr[0][2], str(d))
        plt.plot(arr[1][1], arr[1][2],'.g', markersize=10)
        plt.plot(arr[2][1], arr[2][2],'.g', markersize=10)
    plt.show()
