# common conversion point stacking

import CCP_plottingroutines as CCP_plot
import matplotlib.pyplot as plt
import sys
import matplotlib.pylab as pylab


params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

name = 'CCP_South_Africa'
rffilter = 'jgf1'
conversion = 'prem'
factor = 2.0


CCP = CCP_plot.ccp_volume()
CCP.load_latest(
    name=name,
     filter=rffilter,
     conversion=conversion,
     factor=factor)

print('done loading')




CCP.plot_datacoverage(660,name=name,filter=rffilter, conversion=conversion, factor=factor)
#CCP.plot_crosssection('NS',25,amplify=1.,name=name,filter=rffilter, conversion=conversion, factor=factor,zoom=False, mincoverage=2.)

#CCP.plot_crosssection_any(lon1=-157,lon2=-140,lat1=69,lat2=57,numpoints=42,amplify=0.4,mincoverage=40.,conversion=conversion, zoom=True)#AA
#CCP.plot_crosssection_any(lon1=-159,lon2=-131,lat1=66.5,lat2=59.5,numpoints=42,amplify=0.4,mincoverage=40., conversion=conversion,  zoom=True)#BB
#CCP.plot_crosssection_any(lon1=-160,lon2=-146,lat1=63,lat2=59,numpoints=24,amplify=0.4,mincoverage=40., conversion=conversion, zoom=True)#CC
#CCP.plot_crosssection_any(lon1=-164,lon2=-154,lat1=60.5,lat2=56.5,numpoints=24,amplify=0.4,mincoverage=40.,convesion=conversion, zoom=True)#DD


#CCP.plot_topography(380, 440,name=name,filter=rffilter,conversion=conversion,factor=factor,mincoverage=2., amplitude = False, blobs = False)
#CCP.plot_mtzwidth(filter=rffilter, conversion=conversion, factor=factor, mincoverage=2.)
#CCP.plot_mtzwidth_write(filter=rffilter, conversion=conversion, factor=factor,mincoverage=2.)
#CCP.plot_moveout(d660=False,filter=rffilter, conversion=conversion, factor=factor)

plt.show()
