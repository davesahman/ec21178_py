#!/usr/bin/env python3

import sys
import numpy as np
from numpy import exp, linspace,random, math
from astropy.stats import LombScargle
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
import glob
from hipercam.hlog import Hlog
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import gaussian_fwhm_to_sigma

np.set_printoptions(precision=12)
lcfile     = 'ydiff.log'
Title      = 'YDIFF'
Norm_start = 0.
device     = '1/xs'
period     = True

x,y,e = np.loadtxt(lcfile, unpack=True, usecols=(0,1,2))

plt.figure(figsize=(20,10))
plt.plot(x,y)
plt.show()

if period:
    # x *= 1440.                                 # Change to minutes
    ls    = LombScargle(x,y,e)                 # Create periodogram
    fmax  = 1                               # Set upper frequency (cycles/min) limit 
    nfreq = int(1000*fmax*(x.max()-x.min()))     # Calculate number of frequency steps, oversample x10
    freq  = np.linspace(fmax/nfreq,fmax,nfreq) # Create frequency array
    #freq = np.linspace(0.1,5,100)
    power = ls.power(freq)                      # Calculate periodogram powers            
    fmax  = freq[power==power.max()]        # Calculate peak in periodogram

    # Plot Periodogram
    plt.figure(figsize=(20,10))
    plt.plot(freq, power)
    plt.show()

print("Peak in periodogram at cycles / day. Period of days",1/fmax) 
