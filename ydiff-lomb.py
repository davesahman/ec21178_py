#!/usr/bin/env python3

import sys
import numpy as np
from numpy import exp, linspace,random, math
from astropy.stats import LombScargle
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
from scipy import interpolate
import glob
from hipercam.hlog import Hlog
from hipercam.hlog import Tseries
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import gaussian_fwhm_to_sigma

np.set_printoptions(precision=12)
lcfile     = 'ydiff.log'
Title      = 'YDIFF'
Norm_start = 0.
device     = '1/xs'
period     = True

x,y,e = np.loadtxt(lcfile, unpack=True, usecols=(0,1,2))
'''
plt.figure(figsize=(20,10))
plt.plot(x,y)
plt.show()
'''
# Subtract linear fit

z = np.polyfit(x, y, 2)
print('z',z)
sys.exit()
p = np.poly1d(z)
y -= p(x)

'''
plt.figure(figsize=(20,10))
plt.plot(x,y)
plt.show()
'''

'''
if period:
    # x *= 1440.                                 # Change to minutes
    ls    = LombScargle(x,y,e)                 # Create periodogram
    fmax  = 20                              # Set upper frequency (cycles/min) limit 
    nfreq = int(1000*fmax*(x.max()-x.min()))     # Calculate number of frequency steps, oversample x10
    freq  = np.linspace(fmax/nfreq,fmax,nfreq) # Create frequency array
    #freq = np.linspace(0.1,5,100)
    power = ls.power(freq)                      # Calculate periodogram powers            
    fmax  = freq[power==power.max()]        # Calculate peak in periodogram

    # Plot Periodogram
    
    plt.figure(figsize=(20,10))
    plt.plot(freq, power)
    plt.show()
    '''

# print("Peak in periodogram at cycles / day. Period of days",1/fmax) 


mask = np.zeros_like(x).astype('int')
ts = Tseries(x,y,e,mask)
period = 0.15452941957207128
t1 = 0.196946800588
ts2 = ts.fold(period,t1)
ts2.y = np.absolute(ts2.y)
ts3 = ts2.bin(400,'mean')


plt.figure(figsize=(16,8))
# plt.scatter(ts2.t,ts2.y,s=4)
plt.plot(ts3.t,ts3.y)
plt.plot(ts3.t+1.0,ts3.y)
plt.xlim(-0.5,1.5)
# plt.plot(ts.t,ts.y)
# plt.plot(ts.t,ynew)
# plt.plot(ts3.t,ydiff)

plt.show()
