#!/usr/bin/env python3

import sys
import numpy as np
from numpy import exp, linspace,random, math
from astropy.timeseries import LombScargle
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
x,y,e = np.loadtxt(lcfile, unpack=True, usecols=(0,1,2))

# Subtract linear fit

z = np.polyfit(x, y, 1)
p = np.poly1d(z)
y -= p(x)

# Create hipercam Tseries object

mask = np.zeros_like(x).astype('int')
ts = Tseries(x,y,e,mask)
period = 0.15452941957207128
t1 = 0.196946800588
ts2 = ts.fold(period,t1)   # Fold on the period
# ts2.y = np.absolute(ts2.y) # make all values positive
ts3 = ts2.bin(40,'mean')   # Bin the plot

plt.figure(figsize=(12,6))
plt.plot(ts3.t,ts3.y)
plt.plot(ts3.t+1.0,ts3.y)
plt.xlim(-0.5,1.5)

plt.show()
