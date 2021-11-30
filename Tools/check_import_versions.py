# Script to compare imported packages to a working version.

# Bash command to check unique imports:
# grep "import" ../Processing_Scripts/*py ../Migration_Scripts/*py ../Stacking_Scripts/*py ../Plotting_Scripts/*py | awk -F":" '{print $2}' | sort | uniq


import geographiclib
import cartopy
import cmath
import concurrent.futures
import glob
import math
import matplotlib
import mpl_toolkits
import msgpack
import msgpack_numpy
import numpy
import obspy
import os
import pickle
import scipy
import shapely
import shutil
import subprocess
import sys
import time
import warnings


print()
print('--- Current system package versions ---')
print()
print('python',sys.version)
print()
print('geographiclib',geographiclib.__version__)
print('cartopy',cartopy.__version__)
print('matplotlib',matplotlib.__version__)
print('numpy',numpy.__version__)
print('obspy',obspy.__version__)
print('scipy',scipy.__version__)
print('shapely',shapely.__version__)

############ Packages without version attribute. ##############
# print('cmath',cmath.__version__)
# print('concurrent',concurrent.futures.__version__)
# print('glob',glob.__version__)
# print('math',math.__version__)
# print('mpl_toolkits',mpl_toolkits.__version__)
# print('msgpack',msgpack.__version__)
# print('msgpack_numpy',msgpack_numpy.__version__)
# print('os',os.__version__)
# print('pickle',pickle.__version__)
# print('shutil',shutil.__version__)
# print('subprocess',subprocess.__version__)
# print('sys',sys.__version__)
# print('time',time.__version__)
# print('warnings',warnings.__version__)
print()
print('---------------------------------------')
print()




###################################################################################
#-----------------    Copy of output from above from working version -------------#
###################################################################################

print()
print('--- Checked system package versions ---')
print()

print('python 3.8.8 (default, Apr 13 2021, 12:59:45)')
print('[Clang 10.0.0 ]')

print('geographiclib 1.52')
print('matplotlib 3.3.4')
print('numpy 1.20.1')
print('obspy 1.2.2')
print('scipy 1.6.2')
print('shapely 1.7.1')

print()
print('---------------------------------------')
print()
