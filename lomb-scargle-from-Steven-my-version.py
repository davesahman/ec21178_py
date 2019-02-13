#!/usr/bin/env python3

import numpy as np
from numpy import exp, linspace,random
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
# print("fmax =",fmax)

period_time = 1/fmax
phase=fmax*x
phase=np.mod(phase,1)
t2 = x - (period_time * int(x[0]/period_time)) # subtract integer of first period 
cycle1 = np.floor_divide(t2,period_time)
# print('cycl_1.shape =',cycle1.shape)

# print ('phase  ',phase[:100])
# print("False alarm prob = {:.20f}".format(ls.false_alarm_probability(power.max())))

'''
a = x[55:90]
b = y[55:90]
c = e[55:90]
# print ('phase =',phase[55:90])

# subtract off integer part of MJD
# otherwise floating point errs in minimisation
# tfloor = int(a.min())
# a -= tfloor

# Fit Gaussian to the min value of b
#-------------------------------------------------------
# find min values of b
min_arg = np.argmin(b)
t_min = a[min_arg]

# get gaussian around peak to fit
b -= b.max()
g1 = np.fabs(b)

# Gaussian fit
# Initial guesses
amp = -1 * b[min_arg]
cen = t_min
wid = est_fwhm(a, g1) * gaussian_fwhm_to_sigma
init_vals = [amp ,cen, wid]
best_vals, covar = curve_fit(gaussian, a, g1, p0=init_vals)

# Plot results
# plt.plot(a, gaussian(a, *best_vals), label="fit")
# plt.plot(a, g1, label="data")
# plt.legend()
# plt.show()

mid_ecl = tfloor + best_vals[1]
mid_ecl_err = np.sqrt(covar[1,1])
print("Eclipse Time = ",mid_ecl, "+/-", mid_ecl_err)
'''

# create array with number of cycle
'''
cycle = np.zeros(phase.shape,dtype=int)
cycle[0]=0
# print ('range(len(phase))',range(len(phase)))
for i in range(len(phase)-1):
    cycle[i + 1] = cycle[i]
    if phase[i] > phase[i+1]:
        cycle[i+1] = cycle[i] + 1

# print('cycle.shape =',cycle.shape)
'''
# create list of bad or missed eclipses

bad_ecl = [1]

# loop over each cycle and fit gaussian
ecl = np.empty([cycle1.shape[0],3],dtype=float)

i = 84
# print('cycle1.max =',cycle1.max())
# input("Press enter to continue......")
# while i < (cycle1.max() - 1):
while i < 100:
    i += 1
    if i in cycle1:
      n=i
      print('i before =',i)
      if i == 86 or i ==93:  # this eclipse was not caught in the data
         continue
      print('i after =',i)
      c_range = cycle1 == i
      # print('c_range = ',c_range)
      x_range = x[c_range]
      # print('x_range = ',x_range)
      y_range = y[c_range]
      tfloor = int(x_range.min())
      x_range -= tfloor
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

    # Plot results
   
    plt.plot(x_fwtm, gaussian(x_fwtm, *best_vals), label="fit")
    plt.plot(x_range, yr, label="data")
    plt.legend()
    plt.show()
  
    ecl[i,0] +=i
    ecl[i,1] +=best_vals[1] + tfloor
    ecl[i,2] +=np.sqrt(covar[1,1])
# print('ecl =',ecl[0:n+1])
#  print('ecl[n,1] =',ecl[n,1])
#  print('ecl[1,1] =',ecl[1,1])

ave_period = (ecl[n,1] - ecl[1,1])/(n-1)

# print('ave period = ',ave_period)

plt.plot(ecl[1:n,0],ecl[1:n,1],'ro')
plt.show()