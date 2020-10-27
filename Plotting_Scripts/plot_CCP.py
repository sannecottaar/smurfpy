# common conversion point stacking
#-----------------------------------------------------#
import CCP_plottingroutines as CCP_plot
import matplotlib.pyplot as plt
import os,sys
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
#-----------------------------------------------------#

# Command line help
if len(sys.argv) < 7 or str(sys.argv[1]).lower() == 'help':
    print('\n')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print(sys.argv[0])
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('Description:           Plots CCP stack, including data coverage, discontinuity depths maps and XCs')
    print('Inputs:                name, conversion, filter band, factor, mincoverage, plot_type, plot_params')
    print('Outputs:               Various CCP plots to matplotlib and pdfs)\n')
    print('Usage:                 >> python3 plot_CCP.py name, conversion, filter band, factor, mincoverage, plot_type, plot_params')
    print('Format                 1-3: [str], 4-5: [float], 6: [str], 7-- [int/float]')
    print('Recommended:           >> python3 plot_CCP.py CCP_Global prem jgf1 2.0 2.0 COV 410')
    print('Recommended:           >> python3 plot_CCP.py CCP_Global prem jgf1 2.0 2.0 TOPO 410')
    print('Recommended:           >> python3 plot_CCP.py CCP_Global prem jgf1 2.0 2.0 THICK')
    print('Recommended:           >> python3 plot_CCP.py CCP_Global prem jgf1 2.0 2.0 GETTHICK')
    print('Recommended:           >> python3 plot_CCP.py CCP_Global prem jgf1 2.0 2.0 MOVEOUT 660')
    print('Recommended:           >> python3 plot_CCP.py CCP_Global prem jgf1 2.0 2.0 NS 20 1.5')
    print('Recommended:           >> python3 plot_CCP.py CCP_Global prem jgf1 2.0 2.0 EW 20 1.5')
    print('Recommended:           >> python3 plot_CCP.py CCP_Global prem jgf1 2.0 2.0 XC 18 -38 48 25 1.5 60')
    print('-----------------------------------------------------------------------------------------------------------------------')
    print('\n')
    sys.exit()

name = str(sys.argv[1])   # CCP_volume name
conversion = str(sys.argv[2])   # depth conversion used
rffilter = str(sys.argv[3])  # RF filter
factor = float(sys.argv[4])  #Fresnel Zone smoothing factor
mincoverage = float(sys.argv[5]) # Minimum data coverage (not always required)
plot_type = str(sys.argv[6]) # Determine the plot type:

CCP = CCP_plot.ccp_volume()
CCP.load_latest(
    name=name,
     filter=rffilter,
     conversion=conversion,
     factor=factor)

print('done loading')
# Make sure a CCP_Figures directory is present
savepath='../CCP_Figures/'
if not os.path.exists(savepath):
    os.makedirs(savepath)

# ------------------------ Data coverage ---------------- ####
if plot_type == 'COV':
    if len(sys.argv) != 8:
        print('Try >> python3 plot_CCP.py help')
        sys.exit('Arguments incorrect.... exit')
    plot_depth=int(sys.argv[7])
    print('Plotting data coverage at '+ str(plot_depth) + 'km')
    CCP.plot_datacoverage(plot_depth,name=name, conversion=conversion, filter=rffilter, factor=factor)
    plt.savefig(savepath + name + '_data_coverage_' + str(conversion) + '_' + str(rffilter) + '_' + str(plot_depth) + '.pdf')

# ------------------- MTZ Discontinuity topography ---------------- ####
if plot_type == 'TOPO':
    if len(sys.argv) != 8:
        print('Try >> python3 plot_CCP.py help')
        sys.exit('Arguments incorrect.... exit')
    plot_depth=int(sys.argv[7])
    print('Plotting topography of '+ str(plot_depth) + 'km discontinuity')
    CCP.plot_topography(plot_depth-40, plot_depth+80,name=name, conversion=conversion, filter=rffilter, factor=factor, mincoverage=2., amplitude = False, blobs = False)
    plt.savefig(savepath + name + '_discont_topo_' + str(conversion) + '_' + str(rffilter) + '_' + str(plot_depth) + '.pdf')

# --------------------  MTZ thickness         ---------------- ####
if plot_type == 'THICK':
    if len(sys.argv) != 7:
        print('Try >> python3 plot_CCP.py help')
        sys.exit('Arguments incorrect.... exit')
    print('Plotting MTZ thickness')    
    CCP.plot_mtzwidth(name=name, conversion=conversion, filter=rffilter, factor=factor, mincoverage=2.)
    plt.savefig(savepath + name + '_MTZ_thickness_' + str(conversion) + '_' + str(rffilter) + '.pdf')

# ----------------  Getting MTZ thickness file ---------------- ####
if plot_type == 'GETTHICK':
    if len(sys.argv) != 7:
        print('Try >> python3 plot_CCP.py help')
        sys.exit('Arguments incorrect.... exit')
    print('Getting MTZ thickness')
    # Get MTZ thickness and save to file at ../CCP_volumes/*txt
    CCP.plot_mtzwidth_write(name=name, conversion=conversion, filter=rffilter, factor=factor, mincoverage=2.)

# ----------------  CCP volume moveout plot ---------------- ####
if plot_type == 'MOVEOUT':
    if len(sys.argv) != 8:
        print('Try >> python3 plot_CCP.py help')
        sys.exit('Arguments incorrect.... exit')
    plot_depth=int(sys.argv[7])
    print('Plotting moveout aligned on '+ str(plot_depth) + 'km ')
    if plot_depth == int(660):
        sort_depth = True
    else:
        sort_depth = False
    CCP.plot_moveout(d660=sort_depth,name=name, conversion=conversion, filter=rffilter, factor=factor)
    plt.savefig(savepath + name + '_moveout_plot_aligned_' + str(conversion) + '_' + str(rffilter) + '_' + str(plot_depth) + '.pdf')

# ----------------  CCP volume cross section NS ---------------- ####
if plot_type == 'NS':
    if len(sys.argv) != 9:
        print('Try >> python3 plot_CCP.py help')
        sys.exit('Arguments incorrect.... exit')
    lonorlat_seed=int(sys.argv[7])
    amp=float(sys.argv[8])
    print('Plotting NS XC at '+ str(lonorlat_seed) + ' Deg longitude')        
    CCP.plot_crosssection(str(plot_type),lonorlat_seed,name=name, conversion=conversion, filter=rffilter, factor=factor, mincoverage=2.,amplify=amp,zoom=True)
    plt.savefig(savepath + name + '_XC_' + str(plot_type) + '_at_' + str(lonorlat_seed) + 'Deg_' + str(conversion) + '_' + str(rffilter) + '.pdf')

# ----------------  CCP volume cross section EW ---------------- ####
if plot_type == 'EW':
    if len(sys.argv) != 9:
        print('Try >> python3 plot_CCP.py help')
        sys.exit('Arguments incorrect.... exit')
    lonorlat_seed=int(sys.argv[7])
    amp=float(sys.argv[8])
    print('Plotting EW XC at '+ str(lonorlat_seed) + ' Deg latitude')     
    CCP.plot_crosssection(str(plot_type),lonorlat_seed,name=name, conversion=conversion, filter=rffilter, factor=factor, mincoverage=2.,amplify=amp,zoom=True)
    plt.savefig(savepath + name + '_XC_' + str(plot_type) + '_at_' + str(lonorlat_seed) + 'Deg_' + str(conversion) + '_' + str(rffilter) + '.pdf')

# -----------  CCP volume cross section Any direction   -------- ####
if plot_type == 'XC':
    if len(sys.argv) != 13:
        print('Try >> python3 plot_CCP.py help')
        sys.exit('Arguments incorrect.... exit')
    print('Plotting defined profile:')
    profile = [int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9]), int(sys.argv[10])]
    amp=float(sys.argv[11])
    np=int(sys.argv[12])
    print('   From '+ str(profile[0]) + 'E,'+str(profile[1])+'N to '+str(profile[2])+'E,'+str(profile[3])+'N')
    CCP.plot_crosssection_any(lon1=profile[0],lon2=profile[2],lat1=profile[1],lat2=profile[3],name=name, conversion=conversion, filter=rffilter, factor=factor, mincoverage=2., amplify=amp, numpoints=np, zoom=True)
    plt.savefig(savepath + name + '_XC_at_' +str(profile[1])+'N'+str(profile[0])+'E_to_' + str(profile[3])+'N'+str(profile[2])+'E_' + str(conversion) + '_' + str(rffilter) + '.pdf')
    
plt.show()
