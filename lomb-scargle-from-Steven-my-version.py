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

# Set up gaussian model
def gaussian(x, amp, cen, wid):
    return amp * exp (-(x-cen)**2/(2*wid**2))

def est_fwhm(x, y):
    """
    Estimate FWHM of a Gaussian
    """
    half_max = 0.5*y.max()
    within_halfmax = y > half_max
    x_within_halfmax = x[within_halfmax]
    return x_within_halfmax.max() - x_within_halfmax.min()

def est_fwtm(x, y):
    """
    Estimate x range of a Gaussian - within tenth of max
    """
    fwtm = 0.1*y.max()
    within_max = y > fwtm
    x_fwtm = x[within_max]
    y_fwtm = y[within_max]
    return x_fwtm, y_fwtm

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return a.flat[idx]

MJD        = 58324
lcfile     = '../data/EC21178-5417.txt'
Title      = 'EC21178'
Norm_start = 0.701
device     = '1/xs'
period     = True

x,y,e = np.loadtxt(lcfile, unpack=True, usecols=(0,1,2))

'''
plt.plot(x,y)
plt.show()
'''

#  Fit Lomb-Scargle periodogram to get approx period

if period:
    # x *= 1440.                                 # Change to minutes
    ls    = LombScargle(x,y,e)                 # Create periodogram
    fmax  = 30.                                 # Set upper frequency (cycles/min) limit 
    nfreq = int(1000*fmax*(x.max()-x.min()))     # Calculate number of frequency steps, oversample x10
    freq  = np.linspace(fmax/nfreq,fmax,nfreq) # Create frequency array
    #freq = np.linspace(0.1,5,100)
    power = ls.power(freq)                      # Calculate periodogram powers            
    fmax  = freq[power==power.max()]        # Calculate peak in periodogram
# print("Peak in periodogram at cycles / day. Period of days",1/fmax) 
period_time = 1/fmax
phase = fmax*x
phase = np.mod(phase,1)
tfloor = period_time * int(x[0]/period_time)
# tfloor = int(x[0])
print('tfloor =',tfloor)
x -= tfloor # subtract integer of first period 
cycle1 = np.floor_divide(x,period_time)
cycle_vals = np.unique(cycle1)
# print('cycle_vals',cycle_vals)
# print('cycle_vals.max()',cycle_vals.max())

# loop over each cycle and fit gaussian

ecl = np.zeros([int(cycle_vals.max()),3],dtype=float)

i = 0
while i < (cycle_vals.max()-1):
# while i < 100:
    i += 1
    if i in cycle1:
      # print('processing eclipse no. =',i)
      if i == 86 or i ==93 or i ==149 or i==152:  # these eclipses are bad
         continue
      c_range = cycle1 == i
      x_range = x[c_range]
      y_range = y[c_range]
      #tfloor = int(x_range.min())
      #x_range -= tfloor
      min_arg = np.argmin(y_range)
      t_min = x_range[min_arg]
      y_range -= y_range.max()
      yr = np.fabs(y_range)
      fwtm = 0.4*yr.max()
      within_max = yr > fwtm
      x_fwtm = x_range[within_max]
      y_fwtm = yr[within_max]

    # Gaussian fit
    # Initial guesses
      amp = yr[min_arg]
      cen = t_min
      wid = est_fwhm(x_fwtm, y_fwtm) * gaussian_fwhm_to_sigma
      init_vals = [amp ,cen, wid]
      best_vals, covar = curve_fit(gaussian, x_fwtm, y_fwtm, p0=init_vals)

    # Plot results of gaussian fits
      '''
      plt.plot(x_fwtm, gaussian(x_fwtm, *best_vals), label="fit")
      plt.plot(x_range, yr, label="data")
      plt.legend()
      # plt.show()
      '''
      ecl[i,0] +=i
      ecl[i,1] +=best_vals[1] # + tfloor
      ecl[i,2] +=np.sqrt(covar[1,1])
      print('ecl[i,0] & ecl[i,1]',ecl[i,0],ecl[i,1])
      

    # scale y values to y.min
      '''
      print('y[c_range]',y[c_range])
      print('y.min()',y.min())
      print('y[c_range].min()',y[c_range].min())
      
      y[c_range] -= (y[c_range].min() - y.min()) 
      '''

# remove rows in ecl array with all zero entries
ecl = ecl[~np.all(ecl==0, axis=1)]
s = ecl.shape[0]
ave_period = (ecl[s-1,1] - ecl[0,1])/(ecl[s-1,0]-1)

print('ave period (days)          = ',ave_period)
print('ave peri(od (mins)          = ',1440*ave_period)
print('Lomb Scargle period (mins) = ',1440*period_time)


# plot eclipse times vs time

plt.figure(figsize=(20,10))
plt.plot(ecl[0:s-1,1],ecl[0:s-1,0], 'ro')
plt.xlim(0.1,28.) 
plt.show()


# Create revised phase array

phase1 = x/ave_period
phase1 = np.mod(phase1,1)
'''
plt.scatter(phase1,y, s=1, marker='o')
# plt.plot(phase1,y,'-', lw=0.4 )
plt.show()
'''
# np.savetxt("eclipse_times.txt",ecl)

exit()