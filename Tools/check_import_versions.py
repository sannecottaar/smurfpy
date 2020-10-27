# Script to compare imported packages to a working version.

# Bash command to check unique imports:
# grep "import" ../Processing_Scripts/*py ../Migration_Scripts/*py ../Stacking_Scripts/*py ../Plotting_Scripts/*py | awk -F":" '{print $2}' | sort | uniq


import geographiclib
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
print('python 3.6.9 (default, Oct  8 2020, 12:12:24)')
print('[GCC 8.4.0]')
print()
print('geographiclib 1.49')
print('matplotlib 3.0.3')
print('numpy 1.16.3')
print('obspy 1.1.1')
print('scipy 1.2.1')
print('shapely 1.6.4.post2')

print()
print('---------------------------------------')
print()