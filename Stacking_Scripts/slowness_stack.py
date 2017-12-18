#
import obspy
from obspy.core import trace
from obspy import read

print(obspy.__version__)
import matplotlib.pyplot as plt
import os.path
import time
import glob
import shutil
import numpy as np
#import receiver_function as rf
import subprocess
import scipy
from scipy import interpolate
import pickle
import rosediagram
import gc





##########################################################################################

def compute_fromlist(stalist,filefig,fileout,stacklat,stacklon,stackrad,comp='jgf1',dcorr=None):
    refdepths=[]
    reftimes=[]
    reftakeoff=[]
    refslow=[]
    azimuths=[]
    eqdepths=[]

    # referencing to 65 degrees
    test=['taup_time -mod prem -deg 65 -h 20 -ph P,Pms,P210s,P410s,P660s']
    out=subprocess.check_output(test,shell=True).decode('utf-8')
    t= out.split()
    l=[x for x in range(len(t)) if t[x]=='P']
    Pt= float(t[l[0]+1])
    Pslow=float(t[l[0]+2])
    Pi=float(t[l[0]+3])
    refdepths.append(0.)
    reftimes.append(0.)
    refslow.append(0.)
    reftakeoff.append(np.nan)
    l=[x for x in range(len(t)) if t[x]=='Pms']
    refdepths.append(33.)
    reftimes.append( float(t[l[0]+1]))
    refslow.append(float(t[l[0]+2])-Pslow)
    reftakeoff.append(float(t[l[0]+3]))
    l=[x for x in range(len(t)) if t[x]=='P210s']
    refdepths.append(220.)
    reftimes.append( float(t[l[0]+1])-Pt)
    refslow.append(float(t[l[0]+2])-Pslow)
    reftakeoff.append(float(t[l[0]+3]))
    l=[x for x in range(len(t)) if t[x]=='P410s']
    refdepths.append(400.)
    reftimes.append( float(t[l[0]+1])-Pt)
    refslow.append(float(t[l[0]+2])-Pslow)
    reftakeoff.append(float(t[l[0]+3]))
    l=[x for x in range(len(t)) if t[x]=='P660s']
    refdepths.append(670.)
    reftimes.append( float(t[l[0]+1])-Pt)
    refslow.append(float(t[l[0]+2])-Pslow)
    reftakeoff.append(float(t[l[0]+3]))

    timespace=np.linspace(-10,120,1210)



    slow_poly=np.polyfit(reftimes,refslow,4)
    slowness=np.poly1d(slow_poly)
    depths_poly=np.polyfit(reftimes,refdepths,4)
    depths=np.poly1d(depths_poly)



    f=plt.figure(figsize=(8,12)) 


    slow=np.linspace(-0.7,0.7,200)
    depthsspace=np.linspace(-10,1500,1510)

    c=0
    # Download all data too list

    for i in range(len(stalist)): #range(cat.count()):

        if stalist[i][:3]=='Dat':
            stalist[i]='../'+stalist[i]
        print('stacking', i,stalist[i])

        if os.path.isfile(stalist[i]):

            seis=read(stalist[i],format='PICKLE')

            dt=0.1#seis[0].stats['delta']
            dist=seis[0].stats['dist']
            baz=seis[0].stats['baz']

            RF=trace.Trace()
            # load receiver function
            if comp=='rff2':
                RF.data= np.real(seis[0].rff2['iterativedeconvolution'])
            if comp=='rff1':
                RF.data= np.real(seis[0].rff1['iterativedeconvolution'])
            if comp=='jgf1':
                RF.data= np.real(seis[0].jgf1['iterativedeconvolution'])
            if comp=='jgf2':
                RF.data= np.real(seis[0].jgf2['iterativedeconvolution'])
            time=seis[0].jgf1['time']
            RF.data=RF.data/np.max(np.abs(RF.data))

            indm=np.argmax(np.abs(RF.data[200:260]))+200

            RF.taper(max_percentage=0.05, type='cosine')

            RFtmp=np.zeros((len(RF.data)+int(2.*20./dt),))
            b=int(20/dt-indm+25/dt)
            if np.mean(RF.data[indm-10:indm+10])<0.:
                RF.data=RF.data*-1.
            print(len(RF.data))
            RFtmp[b:b+len(RF.data)]= RF.data/np.max(np.abs(RF.data))

            RFdepth=interpolate.griddata(seis[0].conversions['prem']['depthsfortime'],RF.data,depthsspace)
 
            # build slowness and depth-converted stack
            if (c==0):
                    slowstack=np.zeros((len(slow),len(RF.data))  )   
                    depthstack=np.zeros(len(RFdepth),)
                    timeref=time
                    if dcorr:
                        depthstackcorr=np.zeros(len(depthsspace),)
                        nancount =np.zeros(len(depthsspace),)
            # stack after depth conversion
            depthstack=depthstack+RFdepth

            # stack 3D corrected data
            if dcorr:
                RFdepthcorr=interpolate.griddata(seis[0].conversions[dcorr]['depthsfortime'],RF.data,depthsspace)
                nancount[~np.isnan(RFdepthcorr)]=nancount[~np.isnan(RFdepthcorr)]+1.
                RFdepthcorr[np.isnan(RFdepthcorr)]=0
                depthstackcorr=depthstackcorr+RFdepthcorr

            # stack in slowness domain
            for s in range(len(slow)):
                    timeshift=int(round((65.-dist)*slow[s]/dt))
                    if (len(RFtmp[int(20./dt-timeshift):int(20./dt-timeshift+np.shape(slowstack)[1])])==np.shape(slowstack[1])[0]):
                        slowstack[s,:]=slowstack[s,:]+RFtmp[int(20./dt-timeshift):int(20./dt-timeshift+np.shape(slowstack)[1])]
	   
            #if plot:
            #    plt.subplot(411)
            #    dist=(seis[0].stats['distancedg'])
            #    plt.scatter(time,dist*np.ones(np.shape(time)),s=1,c=RF.data,vmin=-.5,vmax=.5, edgecolors='none')

 
            c=c+1


    # calculate normalized stack
    slowstack=slowstack/float(c)
    slowstacknorm=slowstack.copy()
    int35=int(35./dt)
    int55=int(55./dt)
    slowstacknorm[:,0:int55]=slowstack[:,0:int55]/np.nanmax(np.abs(slowstack[:,0:int55]))
    slowstacknorm[:,35./dt:55./dt]=slowstack[:,35./dt:55./dt]/np.nanmax(np.abs(slowstack[:,35./dt:55./dt]))
    slowstacknorm[:,int55:-1]=slowstack[:,int55:-1]/np.nanmax(np.abs(slowstack[:,int55:-1])) 


    slowstackmax1=np.nanmax(np.abs(slowstack[:,int35:int55]))/np.nanmax(np.abs(slowstack[:,0:int35]))
    slowstackmax2=np.nanmax(np.abs(slowstack[:,int55:-1]))/np.nanmax(np.abs(slowstack[:,0:int35]))
    #normalize depthstack
    x1=np.argmin(np.abs(depthsspace-depths(10.)))
    x2=np.argmin(np.abs(depthsspace-depths(30.)))

    depthstack=depthstack/c
    depthstacknorm=depthstack.copy()   
    depthstacknorm[0:x1]=depthstack[0:x1]/np.nanmax(np.abs(depthstack[0:x1]))
    depthstacknorm[x1:x2]=depthstack[x1:x2]/np.nanmax(np.abs(depthstack[x1:x2]))
    depthstacknorm[x2:-1]=depthstack[x2:-1]/np.nanmax(np.abs(depthstack[x2:-1])) 

    depthstackmax1=np.nanmax(np.abs(depthstack[x1:x2]))/np.nanmax(np.abs(depthstack[0:x1]))
    depthstackmax2=np.nanmax(np.abs(depthstack[x2:]))/np.nanmax(np.abs(depthstack[0:x1]))

    if dcorr:
        depthstackcorr=depthstackcorr/nancount
        depthstackcorrnorm=depthstackcorr.copy()   
        depthstackcorrnorm[0:x1]=depthstackcorr[0:x1]/np.nanmax(np.abs(depthstackcorr[0:x1]))
        depthstackcorrnorm[x1:x2]=depthstackcorr[x1:x2]/np.nanmax(np.abs(depthstackcorr[x1:x2]))
        depthstackcorrnorm[x2:-1]=depthstackcorr[x2:-1]/np.nanmax(np.abs(depthstackcorr[x2:-1])) 

    # take RF out of slowness stack
    xx,yy=np.meshgrid(time,slow)
    RF_timestack=interpolate.griddata((xx.ravel(),yy.ravel()),slowstack.ravel(),(timespace,slowness(timespace)),method='cubic')
    RF_timestacknorm=interpolate.griddata((xx.ravel(),yy.ravel()),slowstacknorm.ravel(),(timespace,slowness(timespace)),method='cubic')
    t1=int(np.argmin(np.abs(timespace-(10.))))
    t2=int(np.argmin(np.abs(timespace-(30.))))
    RF_timestacknorm= RF_timestack.copy()

    RF_timestacknorm[0:t1]=RF_timestacknorm[0:t1]/np.max(np.abs(RF_timestack[0:t1]))
    RF_timestacknorm[t1:t2]=RF_timestacknorm[t1:t2]/np.max(np.abs(RF_timestack[t1:t2]))
    RF_timestacknorm[t2:-1]=RF_timestacknorm[t2:-1]/np.max(np.abs(RF_timestack[t2:-1]))
    timestackmax1=np.max(np.abs(RF_timestack[t1:t2]))/np.max(np.abs(RF_timestack[0:t1]))
    timestackmax2=np.max(np.abs(RF_timestack[t2:-1]))/np.max(np.abs(RF_timestack[0:t1]))
    
    RF_depthstack=interpolate.griddata(depths(timespace),RF_timestack,depthsspace)
    RF_depthstacknorm=interpolate.griddata(depths(timespace),RF_timestacknorm,depthsspace)

    c=0



    # calculate errors to stacks   
    for i in range(len(stalist)): #range(cat.count()):
        print('calculating errors', i,stalist[i])
        if os.path.isfile(stalist[i]):
            seis=read(stalist[i],format='PICKLE')
            if comp=='rff1':
                ref=seis[0].rff1
            elif comp=='rff2':
                ref=seis[0].rff2
            elif comp=='jgf1':
                ref=seis[0].jgf1
            elif comp=='jgf2':
                ref=seis[0].jgf2
            RF.data= np.real(ref['iterativedeconvolution'])
 
            indm=np.argmax(np.abs(RF.data[200:260])) +200
            RF.data=RF.data/np.max(np.abs(RF.data))
            RFdepth=interpolate.griddata(seis[0].conversions['prem']['depthsfortime'],RF.data,depthsspace)
            RFtmp=np.zeros(int((len(RF.data)+2.*20./dt),))

            RF.taper(max_percentage=0.05, type='cosine')
            b=int(20/dt-indm+25/dt)
            if np.mean(RF.data[indm-10:indm+10])<0.:
                RF.data=RF.data*-1.
            RFtmp[b:b+len(RF.data)]= RF.data/np.max(np.abs(RF.data))

            b=int(20/dt-indm+25/dt)
            RFtmp[b:b+len(RF.data)]= RF.data/np.max(np.abs(RF.data))
            RFtmpnorm=RFtmp.copy()
            int35=int(35./dt)
            int55=int(55./dt)
            RFtmpnorm[0:int35]= RFtmpnorm[0:int35]/np.max(RFtmp[0:int35])
            RFtmpnorm[int35:int55]= RFtmpnorm[int35:int55]/np.max(RFtmp[int35:int55])
            RFtmpnorm[int55:-1]= RFtmpnorm[int55:-1]/np.max(RFtmp[int55:-1])



            if c==0:
                sigma_depth=np.zeros(len(depthstack),)
            sigma_depth=sigma_depth+np.square(RFdepth-depthstack)


            if (c==0):
                sigma_slowstack=np.zeros((len(slow),len(RF.data)))   
            for s in range(len(slow)):
                timeshift=int(round((65.-dist)*slow[s]/dt))
                try:
                    sigma_slowstack[s,:]=sigma_slowstack[s,:]+np.square(RFtmp[20/dt-timeshift:20/dt-timeshift+np.shape(sigma_slowstack)[1]]-slowstack[s,:])

                except:
                    pass
            if dcorr:
                if c==0:
                    sigma_depthcorr=np.zeros(len(depthstack),)
                RFdepthcorr=interpolate.griddata(seis[0].conversions[dcorr]['depthsfortime'],RF.data,depthsspace)
                sigma_depthcorr=sigma_depthcorr+np.square(RFdepthcorr-depthstack)
            c=c+1

    #normalize and interpolate along slowness prediction
    sigma_slowstack=np.sqrt(sigma_slowstack/(c*(c-1)))
    sigma_slowstacknorm=sigma_slowstack.copy()
    int35=int(35./dt)
    int55=int(55./dt)
    int110=int(110./dt)
    sigma_slowstacknorm[:,0:int35]=sigma_slowstacknorm[:,0:int35]/np.max(np.abs(slowstack[:,0:int35]))
    sigma_slowstacknorm[:,int35:int55]=sigma_slowstacknorm[:,int35:int55]/np.max(np.abs(slowstack[:,int35:int55]))
    sigma_slowstacknorm[:,int55:int110]=sigma_slowstacknorm[:,int55:int110]/np.max(np.abs(slowstack[:,int55:int110]))
    sigma_slowstacknorm[:,int110:-1]=sigma_slowstacknorm[:,int110:-1]/np.max(np.abs(slowstack[:,int110:-1])) 
    sigma_timestack=interpolate.griddata((xx.ravel(),yy.ravel()),sigma_slowstack.ravel(),(timespace,slowness(timespace)),method='cubic')
    #sigma_timestack_norm=interpolate.griddata((xx.ravel(),yy.ravel()),sigma_slowstacknorm.ravel(),(timespace,slowness(timespace)),method='cubic')
    sigma_timestack_norm=sigma_timestack.copy()
    sigma_timestack_norm[0:t1]=sigma_timestack_norm[0:t1]/np.max(np.abs(RF_timestack[0:t1]))
    sigma_timestack_norm[t1:t2]=sigma_timestack_norm[t1:t2]/np.max(np.abs(RF_timestack[t1:t2]))
    sigma_timestack_norm[t2:-1]=sigma_timestack_norm[t2:-1]/np.max(np.abs(RF_timestack[t2:-1]))
  

    if dcorr is not None:
        sigma_depthcorr=np.sqrt(sigma_depthcorr/(c*(c-1)))
        sigma_depthcorr_norm=sigma_depthcorr.copy()
        sigma_depthcorr_norm[0:x1]=sigma_depthcorr_norm[0:x1]/np.max(np.abs(depthstackcorr[0:x1]))
        sigma_depthcorr_norm[x1:x2]=sigma_depthcorr_norm[x1:x2]/np.max(np.abs(depthstackcorr[x1:x2]))
        sigma_depthcorr_norm[x2:-1]=sigma_depthcorr_norm[x2:-1]/np.max(np.abs(depthstackcorr[x2:-1]))


    sigma_depth=np.sqrt(sigma_depth/(c*(c-1)))

    sigma_depth_norm=sigma_depth.copy()
    sigma_depth_norm[0:x1]=sigma_depth_norm[0:x1]/np.max(np.abs(depthstack[0:x1]))
    sigma_depth_norm[x1:x2]=sigma_depth_norm[x1:x2]/np.max(np.abs(depthstack[x1:x2]))
    sigma_depth_norm[x2:-1]=sigma_depth_norm[x2:-1]/np.max(np.abs(depthstack[x2:-1]))


    '''
    plt.subplot(411)
    #plt.xlabel('time(s) from 40s before P wave arrival')
    plt.ylabel('distance(dg)')
    plt.title(str(c)+ ' ' + fileout)
    plt.ylim([35,95])
    plt.xlim([-25,120.])

    # plot backazimuth distibution
    ax=f.add_axes([0.66, 0.65, 0.3, 0.14],projection='polar')
    z = np.cos(np.radians(azimuths))
    coll,counts = rosediagram.rose(azimuths, color_by='count', bins=18)
    plt.xticks(np.radians(range(0, 360, 45)), 
               ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
    try:
        plt.rgrids(range(int(np.mean(counts)),int(max(counts)),int(np.mean(counts))), angle=290)
    except:
        print 'rose diagram failed'
    '''

    plt.subplot(411)
    contour_levels = np.arange(-1.,1.01,.1)
    cs3=plt.contourf(timeref,slow,slowstacknorm,contour_levels,rasterized=True, cmap='seismic',extend='both')#,extend='both')
    #cs3.cmap.set_over(color=(0.5,0,0,1.))
    #cs3.cmap.set_under(color=(0.0,0,0.5,1.))
    plt.plot([10.,10.],[-0.7,.7],'k')
    plt.plot([30.,30.],[-0.7,.7],'k')
    #plt.plot([85.,85.],[-0.7,.3],'k')
    plt.plot(reftimes[0]-reftimes[0],refslow[0]-refslow[0],'k',marker='+',markersize=10)
    plt.plot(reftimes[3]-reftimes[0],refslow[3]-refslow[0],'k',marker='+',markersize=10)
    plt.plot(reftimes[4]-reftimes[0],refslow[4]-refslow[0],'k',marker='+',markersize=10)
    plt.plot(timespace,slowness(timespace),'k')

    #plt.text(18.,0.2,str(round(slowstackmax1,2)))
    #plt.text(60.,0.2,str(round(slowstackmax2,2)))
    #plt.text(110.,0.2,str(round(slowstackmax3,2)))
    plt.ylim([-0.7,.7])
    plt.xlim([-10,120.])
    #plt.colorbar()
    plt.xlabel('time (s) w.r.t. P')
    plt.ylabel('differential slowness (s/deg)')

    plt.subplot(412)
    contour_levels = np.arange(-1.,1.01,.1)
    cs3=plt.contourf(timeref,slow,2.*slowstacknorm,contour_levels,extend='both',rasterized=True, cmap='seismic')
    cs3.cmap.set_over(color=(0.5,0,0,1.))
    cs3.cmap.set_under(color=(0.0,0,0.5,1.))
    plt.plot([10.,10.],[-0.7,.7],'k')
    plt.plot([30.,30.],[-0.7,.7],'k')
    #plt.plot([85.,85.],[-0.7,.3],'k')
    plt.plot(reftimes[0]-reftimes[0],refslow[0]-refslow[0],'k',marker='+',markersize=10)
    plt.plot(reftimes[3]-reftimes[0],refslow[3]-refslow[0],'k',marker='+',markersize=10)
    plt.plot(reftimes[4]-reftimes[0],refslow[4]-refslow[0],'k',marker='+',markersize=10)
    plt.plot(timespace,slowness(timespace),'k')

    #plt.text(18.,0.2,str(round(slowstackmax1,2)))
    #plt.text(60.,0.2,str(round(slowstackmax2,2)))
    #plt.text(110.,0.2,str(round(slowstackmax3,2)))
    plt.ylim([-0.7,.7])
    plt.xlim([-10,120.])
    #plt.colorbar()
    plt.xlabel('time (s) w.r.t. P')
    plt.ylabel('differential slowness (s/deg)')

    plt.subplot(413)
    plt.plot(timespace,RF_timestack,'g',label='normalized')
    plt.plot(timespace,RF_timestacknorm,'b',label='normalized in sections')

    plt.fill_between((timespace),RF_timestacknorm-2.*sigma_timestack_norm,RF_timestacknorm+2.*sigma_timestack_norm,color='b', alpha=0.3)

    plt.xlim([-10,120.])
    plt.ylim([-1.2,1.2])
    plt.plot([10.,10.],[-1.2,1.2],'k')
    plt.plot([30.,30.],[-1.2,1.2],'k')
    plt.plot([-10,120],[0., 0.],'--k')
    plt.plot(reftimes[0]-reftimes[0],[1.0],'k',marker='+',markersize=10)
    plt.plot(reftimes[3]-reftimes[0],[1.0],'k',marker='+',markersize=10)
    plt.plot(reftimes[4]-reftimes[0],[1.0],'k',marker='+',markersize=10)
    plt.xlabel('time (s) w.r.t. P')
    plt.legend(loc='lower right',prop={'size':8})
    # convert time to depth
    plt.subplot(414)
    plt.plot(depths(timespace),RF_timestacknorm,'b', label= 'from slowness stack')
    plt.fill_between(depths(timespace),RF_timestacknorm-2.*sigma_timestack_norm,RF_timestacknorm+2.*sigma_timestack_norm,color='b', alpha=0.3)
    plt.plot(depthsspace,depthstacknorm,'r', label= 'stack of depth-converted RFs')

    plt.fill_between(depthsspace,depthstacknorm-2.*sigma_depth_norm,depthstacknorm+2.*sigma_depth_norm,color='r', alpha=0.3)
    if dcorr:
            plt.plot(depthsspace,depthstackcorrnorm,'g', label= 'stack of 3D corr depth-converted RFs')
            plt.fill_between(depthsspace,depthstackcorrnorm-2.*sigma_depthcorr_norm,depthstackcorrnorm+2.*sigma_depthcorr_norm,color='g', alpha=0.3)
    plt.plot([depths(10.),depths(10.)],[-1.4,1.4],'k')
    plt.plot([depths(30.),depths(30.)],[-1.4,1.4],'k')
    plt.plot([-10,1500],[0., 0.],'--k')
    plt.xlim([-10,1500.])
    plt.ylim([-1.4,1.4])
    plt.xlabel('depth (km)')
    plt.legend(loc='lower right',prop={'size':8})

    #plt.show()

    plt.savefig(filefig)


    # save stack
    stack=dict()

    stack['depthstack']=depthstack
    stack['depthstack_norm']=depthstacknorm
    stack['depthstack_axis']=depthsspace
    stack['depthstack_sigma']=sigma_depth
    stack['depthstack_sigma_norm']=sigma_depth_norm
    if dcorr:
        stack['depthstack_'+dcorr]=depthstackcorr
        stack['depthstack_norm_'+dcorr]=depthstackcorrnorm
        stack['depthstack_axis_'+dcorr]=depthsspace
        stack['depthstack_sigma_'+dcorr]=sigma_depthcorr
        stack['depthstack_sigma_norm_'+dcorr]=sigma_depthcorr_norm

    stack['slowstack_time']=RF_timestack
    stack['slowstack_depth']=RF_depthstack
    stack['slowstack_time_norm']=RF_timestacknorm
    stack['slowstack_depth_norm']=RF_depthstacknorm
    stack['slowstack']=slowstack
    stack['slowstack_norm']=slowstacknorm
    stack['slowstack_timeaxis']=time
    stack['slowstack_slowaxis']=slow
    stack['slowstack_sigma']=sigma_slowstack
    stack['slowstack_sigma_norm']=sigma_slowstacknorm
    stack['slowstack_time_sigma']=sigma_timestack
    stack['slowstack_time_sigma_norm']=sigma_timestack_norm    
    #stack['latitude']=seis[0].stats['latitude']
    #stack['longitude']=seis[0].stats['longitude']
    stack['figure_filename']=filefig
    stack['latitude']=stacklat
    stack['longitude']=stacklon
    stack['radius']=stackrad
    stack['filter']=seis[0].jgf1['filter']
    stack['filterconst']=seis[0].jgf1['filterconst']
    #stack['fitmin']=fitmin
    #stack['ampmax']=ampmax
    stack['stacked_num']=c
    stack['filelist']=stalist
    # write out
    with open(fileout,'wb') as handle:
        pickle.dump(stack,handle)
    plt.show()
 



##########################################################################################




stalist = glob.glob('../Data/*/*PICKLE')

compute_fromlist(stalist, 'test.pdf', 'test.PICKLE',0,0,0)


