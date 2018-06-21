 #

import obspy
from obspy import read
from obspy.core import Stream
#imports
from obspy.core import event

from obspy import UTCDateTime
import matplotlib.pyplot as plt
import os.path
import time
import glob
import shutil
import numpy as np
import scipy
#from scipy import interpolate
#from scipy import stats
import cmath
#import spectrum as sp
#from spectrum.tools import nextpow2
from pylab import semilogy
import warnings
import time






def water_level_decon(comp2, comp1, wl=0.001,dt=0.025,filt='cosine',fmax=0.5,timeshift=25.):
    '''
Water-level deconvolution
After Jennifer Andrew's thesis
Deconvolves comp2 by comp1
wl = the water-level, usually ranges from 1e-5 to 1e-1
a  = width of the gaussian filter
dt = 1/(sampling rate) seconds
filt = 'gauss' or 'cosine'
fmax = gaussian filter width or cut-off frequency for cosine filter
timeshift = time beginning of seismogram and the main arrival in comp1
    ''' 

    print ('Performing water-level deconvolution, with a water-level constant of ', str(wl), ' and a ',filt, ' filter with an fmax of ', str(fmax), ' and the resulting RF is time shifted by ', str(timeshift),' seconds.')


    # pad the reference component

    padded = np.zeros_like(comp2)

    padded[0:comp1.size] = comp1

 
    # convert to frequency domain

    c2freq = np.fft.fft(comp2)

    c1freq = np.fft.fft(padded)

    frq    = np.fft.fftfreq(len(comp2),dt)

 
 
    # construct the denominator, fill in the troughs to stabalize the deconvolution

    wlfreq = c1freq * c1freq.conjugate()
    wlfreqold=wlfreq
    if wl!='search':
        wlfreq[wlfreq.real < wl*wlfreq.real.max()] = complex(wl*wlfreq.real.max(),0.) #water-level


    # Filters
    if filt=='gauss':
        filterf = np.exp( -2.*np.pi*frq*frq/(4.*fmax*fmax) )   #gaussfilter, fmax is the width of the filter and not the maximum frequency present...

    else:
        if filt=='cosine':
            filterf = np.square( np.cos( np.pi*frq/(2.*fmax) ) ) # fmax is the cut-off frequency
            filterf[np.abs(frq)>fmax] = 0.
        else:
            exit('Choose Gauss or Cosine filter for your water-level deconvolution')
            

 

 
    if wl=='search':
        alphas=np.logspace(-4,0,20)
        misfit=np.empty(alphas.shape)
        magntd=np.empty(alphas.shape)
        phaseshift = [cmath.exp(1j*(frq[x]*2.*np.pi*-(timeshift))) for x in range(len(frq))] 
        # find best damping coefficient
        for a in range(len(alphas)):
            wlfreq[wlfreq.real < alphas[a]*wlfreq.real.max()] = complex(alphas[a]*wlfreq.real.max(),0.) #water-level
            RF=np.fft.ifft(c2freq*c1freq.conjugate()/(wlfreq)*filterf*phaseshift,len(comp2))
            comp2c=np.convolve(comp1,RF.real,'same')
            misfit[a]=np.sum(np.square(comp2c-comp2))
            magntd[a]=np.sum(np.abs(RF))
        testalpha=np.square(magntd/np.std(magntd))+np.square(misfit/np.std(misfit))   
        print('best damping factor is ', alphas[testalpha.argmin()])
        wlfreq[wlfreq.real < alphas[testalpha.argmin()]*wlfreq.real.max()] = complex(alphas[testalpha.argmin()]*wlfreq.real.max(),0.) #water-level
        RF=np.fft.ifft(c2freq*c1freq.conjugate()/(wlfreq)*filterf*phaseshift,len(comp2))
    else:
        # Perform deconvolution and filter
        newfreq = c2freq * c1freq.conjugate() / (wlfreq)
        newfreq = newfreq * filterf 
        # Phase shift to recover time
        phaseshift = [cmath.exp(1j* (-frq[x]*2.*np.pi*timeshift) ) for x in range(len(frq))]
        newfreq    = newfreq * phaseshift
        # Convert back to time
        RF  = np.fft.ifft((newfreq),len(comp2))

        # reconvolve with vertical component
        conv=np.real(np.convolve(RF,padded,'full'))
        conv=conv[timeshift/dt:timeshift/dt+len(comp2)]
        # Calculate fit
        comp2f=np.real(np.fft.ifft(filterf*np.fft.fft(comp2)))
        # calculate residual and fit
        residual=comp2f-conv
        fit= 100.* (1.- sum(residual*residual)/sum(comp2f*comp2f))
 
    return RF.real, fit


def multitaper(comp2,comp1,p=2.5, k=4, dt=0.04,filt='cosine',fmax=0.5, timeshift=25.):
    '''
 After Park&Levin 2000, only useful for short receiver functions, (length of taper window)/(2*p)
 Deconvolves comp2 by comp1
 p = time half bandwidth paramater
 k = number of Slepian tapers, k<2p
 filt = 'gauss' or 'cosine'
 fmax = gaussian filter width or cut-off frequency for cosine filter
 timeshift = time beginning of seismogram and the main arrival in comp1

 Returns
 Receiver function in the time domain
 Coherence between the two components in the frequency domain
 Variance of the receiver function in the frequency domain
    '''

    print('Performing multi-taper deconvolution, with a taper bandwidth of ', str(p), ',', str(k),'Slepian tapers,  and a ',filt, ' filter with an fmax of ', str(fmax), ' and the resulting RF is time shifted by ', str(timeshift),' seconds.')
    warnings.warn('Very clumsily stabalized the deconvolution, should be done with a noise estimate or a proporly search of the stabalization factor, like is done in the extended multi-taper function')


    # pad the reference component

    padded_comp1 = np.zeros_like(comp2)

    padded_comp1[:comp1.size] = comp1


    # Compute Slepian data taper

    [w,eigens]=sp.dpss(len(comp2),p,k)  


    # Convert to multi-taper frequency domain

    NFFT = max(256, 2**nextpow2(len(comp2)))

    frq  = np.fft.fftfreq(NFFT,dt)

    Sk2  = np.fft.fft(np.multiply(w.transpose(), comp2.real), NFFT)

    Sk1  = np.fft.fft(np.multiply(w.transpose(), padded_comp1.real), NFFT)   


    # Compute coherence measure in frequency domain (see Park & Levin 2000)

    coherence   =  np.sum(Sk1.conjugate()*Sk2,axis=0)/np.sqrt(np.sum(Sk2.conjugate()*Sk2,axis=0)*np.sum(Sk1.conjugate()*Sk1,axis=0))



    # Filters

    if filt=='gauss':
        filterf = np.exp( -frq*frq/(4.*fmax*fmax) )   #gaussfilter, fmax is the width of the filter and not the maximum frequency present...

    else:

        if filt=='cosine':
            filterf = np.square( np.cos( np.pi*frq/(2.*fmax) ) ) # fmax is the cut-off frequency
            filterf[np.abs(frq)>fmax] = 0.
        else:

            exit('Choose Gauss or Cosine filter for your water-level deconvolution')  

   # compute RF(f)

    denom  =  np.sum(Sk1.conjugate()*Sk1,axis=0)

    RF     =  np.sum(Sk1.conjugate()*Sk2,axis=0)/(denom+0.001*np.max(denom)) # This is stabalized by adding a factor, this could be done better with a noise estimate or finding the best factor, like is done in the extended multi-taper method    

    RF     =  RF*filterf

    # Compute variance in the frequency domain

    variance = (1.-coherence**2.)/((k-1.)*coherence**2) * RF*RF.conjugate()

    # Phase-shift for a time shift in the time domain

    phaseshift = [cmath.exp(1j*(frq[x]*2.*np.pi*-(timeshift))) for x in range(len(frq))]  
    RF         = RF * phaseshift

 
 
    # return to time domain
    RF=np.fft.ifft(RF,len(comp2))



    return RF.real, coherence, variance





def extended_multitaper_noise(comp2,comp1,compnoise,win=50., overlap=0.9, p=2.5, k=4, dt=0.04,filt='cosine',fmax=0.5, timeshift=25.):
    '''
 After Park&Levin 2000 and Helffrich 2006, Shibutani et al. 2008 and and Ved's matlab code (Lekic & Fischer, 2013).
 Different from Ved's code, by summing all the windows for each taper in the frequency domain
 Deconvolves comp2 by comp1
 p = time half bandwidth paramater
 k = number of Slepian tapers, k<2p
 filt = 'gauss' or 'cosine'
 fmax = gaussian filter width or cut-off frequency for cosine filter
 timeshift = time beginning of seismogram and the main arrival in comp1

Choice of p and k can be tested with sum_multi_tapers to gain a flat amplitude range. Don't trust the last 50s of the receiver function. Good choices are (with increasing computational expense)
p=4, k=4, overl=0.75
p=4, k=3, overl=0.9
p=4, k=4, overl=0.9
p=4, k=7, overl=0.95

 Returns
 Receiver function in the time domain

    '''

    print('Performing extended multi-taper deconvolution, with a window width of ', str(win),' and window overlap of ', str(overlap) ,' and  taper bandwidth of ', str(p), ',', str(k),'Slepian tapers,  and a ',filt, ' filter with an fmax of ', str(fmax), ' and the resulting RF is time shifted by ', str(timeshift),' seconds.')




    # pad the reference component
    padded_comp1 = np.zeros_like(comp2)
    padded_comp1[0:comp1.size] = comp1
    padded_compnoise = np.zeros_like(comp2)
    padded_compnoise[0:compnoise.size] = compnoise


    # Compute Slepian data tapers for time windows of 50s with 90% overlap
    lentaper = int(50./dt)
    [Sltapers,eigens] = sp.dpss(lentaper,p,k)   # multi-tapers

    # Define windows t0 use
    starts = np.arange(0,len(comp2)-lentaper,round((1-overlap)*lentaper))

    # Frequency domain
    NFFT = len(comp2)
    frq  = np.fft.fftfreq(NFFT,dt)
    


    # Filters
    if filt=='gauss':
        filterf = np.exp( -frq*frq/(4.*fmax*fmax) )   #gaussfilter, fmax is the width of the filter and not the maximum frequency present...
    else:
        if filt=='cosine':
            filterf = np.square( np.cos( np.pi*frq/(2.*fmax) ) ) # fmax is the cut-off frequency
            filterf[np.abs(frq)>fmax] = 0.
        else:
            exit('Choose Gauss or Cosine filter for your water-level deconvolution') 

    # Sum all time windows for each taper
    denom = np.zeros((NFFT,),dtype=complex)   
    numer = np.zeros((NFFT,),dtype=complex)
    noise = np.zeros((NFFT,),dtype=complex)
    for l in range(k): # Loop over tapers
        # Initialize taper estimates
        SK1  = np.zeros((NFFT,),dtype=complex)
        SK2  = np.zeros((NFFT,),dtype=complex)
        SKN  = np.zeros((NFFT,),dtype=complex)
        SK1 = SK1+(np.fft.fft(np.multiply(Sltapers[:,l], padded_comp1[0:lentaper]), NFFT))
        for n in range(len(starts)): # Loop over windows
            phaseshift = np.array([cmath.exp(1j*(frq[x]*2.*np.pi*-1*starts[n]*dt)) for x in range(len(frq))]  )
            #SK1 = SK1+(np.fft.fft(np.multiply(Sltapers[:,l], padded_comp1[starts[n]:starts[n]+lentaper]), NFFT))*phaseshift
            SK2 = SK2+(np.fft.fft(np.multiply(Sltapers[:,l], comp2[starts[n]:starts[n]+lentaper]), NFFT))*phaseshift
            SKN = SKN+(np.fft.fft(np.multiply(Sltapers[:,l], padded_compnoise[starts[n]:starts[n]+lentaper]), NFFT))*phaseshift
        denom=denom+SK1*SK1.conjugate()   
        numer=numer+SK2*SK1.conjugate()
        noise=noise+SKN*SK1.conjugate()
 

  

    denom[denom.real < noise.real] = noise.real #water-level, should be replaced
    phaseshift = [cmath.exp(1j*(frq[x]*2.*np.pi*-(timeshift))) for x in range(len(frq))] 
    # filter and convert the time domain
    RF=np.fft.ifft(numer/(denom)*filterf*phaseshift,len(comp2))
 

    return RF.real





def extended_multitaper_vs2(comp2,comp1,win=50., overlap=0.9, p=2.5, k=4, dt=0.04,filt='cosine',fmax=0.5, timeshift=25.):
    '''
 After Park&Levin 2000 and Helffrich 2006, Shibutani et al. 2008, and and Ved's matlab code (Lekic & Fischer, 2013)
 This version isn't working properly. The setup also requires the main arrival in comp1 to fall in the center of the first window. 
 Deconvolves comp2 by comp1
 p = time half bandwidth paramater
 k = number of Slepian tapers, k=2p-1
 filt = 'gauss' or 'cosine'
 fmax = gaussian filter width or cut-off frequency for cosine filter
 timeshift = time beginning of seismogram and the main arrival in comp1

Choice of p and k can be tested with sum_multi_tapers to gain a flat amplitude range. Don't trust the last 50s of the receiver function.

 Returns
 Receiver function in the time domain

    '''

    print('Performing extended multi-taper deconvolution, with a window width of ', str(win),' and window overlap of ', str(overlap) ,' and  taper bandwidth of ', str(p), ',', str(k),'Slepian tapers,  and a ',filt, ' filter with an fmax of ', str(fmax), ' and the resulting RF is time shifted by ', str(timeshift),' seconds.')
    warnings.warn('This version of the extenden multi-taper is not working properly, use version 1')


    # Compute Slepian data tapers for time windows of 50s with 90% overlap
    lentaper = int(50./dt)
    [Sltapers,eigens] = sp.dpss(lentaper,p,k)   # multi-tapers

    # Define windows t0 use
    starts = np.arange(0,len(comp2)-lentaper,round((1-overlap)*lentaper))

    # Frequency domain
    NFFT = lentaper
    frq  = np.fft.fftfreq(NFFT,dt)

    # Initialize taper estimates
    SK1  = np.zeros((NFFT,k),dtype=complex)
    SK2  = np.zeros((NFFT,k,len(starts)),dtype=complex)
 
    # Construct denomonimator
    denom = np.zeros((1,NFFT),dtype=complex)   
    for l in range(k): # Loop over tapers
        SK1[:,l] = np.fft.fft( np.multiply(Sltapers[:,l], comp1[0:lentaper]), NFFT) 
        denom = denom + SK1[:,l]*SK1[:,l].conjugate()
    denom[denom.real < .01*denom.real.max()] = .01*denom.real.max() #water-level

    # Construct numerator
    numer = np.zeros((NFFT,len(starts)),dtype=complex)
    for n in range(len(starts)): # Loop over windows
        for l in range(k): # Loop over tapers
            SK2[:,l,n] = (np.fft.fft(np.multiply(Sltapers[:,l], comp2[starts[n]:starts[n]+lentaper]), NFFT)) 
            numer[:,n]=numer[:,n]+SK2[:,l,n]*SK1[:,l].conjugate()

 
    # Filters
    if filt=='gauss':
        filterf = np.exp( -frq*frq/(4.*fmax*fmax) )   #gaussfilter, fmax is the width of the filter and not the maximum frequency present...
    else:
        if filt=='cosine':
            filterf = np.square( np.cos( np.pi*frq/(2.*fmax) ) ) # fmax is the cut-off frequency
            filterf[np.abs(frq)>fmax] = 0.
        else:
            exit('Choose Gauss or Cosine filter for your water-level deconvolution')  

    # Combine all windows in the time domain
    RFwin    = np.zeros((lentaper,len(starts)),dtype=complex)
    RFwinpad = np.empty((len(comp2),len(starts)),dtype=complex)
    RFwinpad.fill(np.nan)
    RF       = np.zeros((len(comp2),1),dtype=complex)
    alphas   = np.logspace(-2,2,20)*np.var(comp2[round(comp2.size/4):3*round(comp2.size/4)])
    for n in range(len(starts)):
        RFtmp      = numer[:,n]/(denom) * filterf
        RFwin[:,n] = np.fft.fftshift(np.fft.ifft(RFtmp.real,lentaper)) 
        RFwinpad[starts[n]:starts[n]+lentaper,n]  =  RFwin[:,n]
    #RF = scipy.stats.nanmean(RFwinpad,axis=1) *1./(1.-overlap) # average all windows
    RF = np.nansum(RFwinpad,axis=1)
    RF  = np.nan_to_num(RF)
 
    return RF



def    sum_multi_tapers(win=2000,lentaper=400,overlap=0.75,p=4,k=7):
    #  Compute Slepian data taper window to test for flat level of amplitudes
        [Sltapers,eigens] = sp.dpss(lentaper,p,k)   # multi-tapers

        # Define windows t0 use
        starts = np.arange(0,win-lentaper,round((1-overlap)*lentaper))
        print(starts)
        sum_multitapers=np.zeros((win,))
        num=np.zeros((win,))
        for n in range(len(starts)):
             for j in range(k):
 
                 print(sum_multitapers[starts[n]:starts[n]+lentaper].shape)
                 sum_multitapers[starts[n]:starts[n]+lentaper]=sum_multitapers[starts[n]:starts[n]+lentaper]+[Sltapers[:,j]]
                 num[starts[n]:starts[n]+lentaper]=num[starts[n]:starts[n]+lentaper]+1.0
        return sum_multitapers
    



def iterative_deconvolution(comp2,comp1, maxbumps=200,dt=0.04,filt='cosine',fmax=0.5, timeshift=25.):
    '''
Ligorria & Ammon 1999, Kikuchi & Kanamori 1982
After Chuck Ammons code 1998- Jenny Andrews 2006
    '''
    print('Performing iterative deconvolution, with a maximum of ',maxbumps,'  and a ',filt, ' filter with an fmax/gausswidth of ', str(fmax), '. The components  should be cut ', timeshift,' seconds before the main arrival.')
 
    # Frequency domain
    NFFT = len(comp2) 
    frq  = np.fft.fftfreq(NFFT,dt)
    
    # Filters
    if filt=='gauss':
        filterf1 = np.exp( -frq*frq/(4.*fmax*fmax) )   #gaussfilter, fmax is the width of the filter and not the maximum frequency present...
        filterf = np.exp( -np.pi*np.pi*frq*frq/(fmax*fmax) )    #Factor of 4 as written in Jenny's script, fortran fft uses angular frequency? 
    else:
        if filt=='cosine':
            filterf = np.square( np.cos( np.pi*frq/(2.*fmax) ) ) # fmax is the cut-off frequency
            filterf[np.abs(frq)>fmax] = 0.
        else:
            exit('Choose Gauss or Cosine filter for your iterative deconvolution')  
    #plt.plot(frq,filterf,label = '1')
    #plt.plot(frq,filterf1,label = '2')
    #plt.plot(frq,filterf2,label = '3')
    #plt.plot([0.5,0.5],[0,1],'k')
    #plt.plot([-5,5],[0.1,.1],'k')
    #plt.plot(frq,np.fft.fft(comp2)/np.max(np.fft.fft(comp2)))
    #plt.plot(frq,np.fft.fft(comp1)/np.max(np.fft.fft(comp1)))
    #plt.legend()
 
    #stophere
    #  phaseshift 
    phaseshift = [cmath.exp(1j*(frq[x]*2.*np.pi*-(timeshift))) for x in range(len(frq))] 

    #Filter components
    comp2f=np.real(np.fft.ifft(filterf*np.fft.fft(comp2)))
    comp1f=np.real(np.fft.ifft(filterf*np.fft.fft(comp1)))
 

    # computer power in numerator
    power=sum(comp1f*comp1f)
    #initialize values
    maxind=np.empty([maxbumps,1],dtype=int)
    maxvalue=np.empty([maxbumps,1],dtype=float)
    residual=comp2f # initial residual is horizontal component
    fit_old=0 # initial fit = 0%

    for peak in range(maxbumps):
    # correlate signals and find peak
        corr=(np.correlate(residual,comp1f,'full'))
        corr=corr[len(comp2)-int(timeshift/dt):-int(timeshift/dt)]# set maximum Xcorr to be at the predicted main arrival
        maxind[peak]=np.argmax(np.absolute(corr))# peak index
        maxvalue[peak]=corr[maxind[peak]]/(power)# peak amplitude


        # build deconvolution result
        decon=np.zeros_like(comp2)
        for ind in range(peak+1):
            decon[maxind[ind]]=decon[maxind[ind]]+maxvalue[ind]

        decon=np.real(np.fft.ifft(filterf*np.fft.fft(decon)))    #filter deconvolution result


        # reconvolve with vertical component
        conv=np.real(np.convolve(decon,comp1,'full'))
        conv=conv[int(timeshift/dt):int(timeshift/dt)+len(comp2)]

        # calculate residual and fit
        residual=comp2f-conv
        fit= 100.* (1.- sum(residual*residual)/sum(comp2f*comp2f))
        if (fit-fit_old)<0.01:
            break
        fit_old=fit

    print('Stopped at', peak, ' when fitting ', fit, ' %')

    return decon, fit
