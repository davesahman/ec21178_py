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



MJD        = 58324
lcfile     = '../data/EC21178-5417.txt'
Title      = 'EC21178'
Norm_start = 0.701
device     = '1/xs'
period     = True

x,y,e = np.loadtxt(lcfile, unpack=True, usecols=(0,1,2))

# x -= MJD

# if len(x[x<Norm_start])>0:
#    e /= y[x<Norm_start].mean()
#    y /= y[x<Norm_start].mean()

# plt.plot(x,y)
# plt.show()

if period:
    # x *= 1440.                                 # Change to minutes
    ls    = LombScargle(x,y,e)                 # Create periodogram
    fmax  = 30.                                 # Set upper frequency (cycles/min) limit 
    nfreq = int(1000*fmax*(x.max()-x.min()))     # Calculate number of frequency steps, oversample x10
    freq  = np.linspace(fmax/nfreq,fmax,nfreq) # Create frequency array
    #freq = np.linspace(0.1,5,100)
    power = ls.power(freq)                      # Calculate periodogram powers            
    fmax  = freq[power==power.max()]        # Calculate peak in periodogram
    print("Peak in periodogram at cycles / min. Period of minutes",1440*1/fmax) 
    print("fmax =",fmax)

phase=fmax*x
# print (phase[:10])
phase=np.mod(phase,1)
# print (phase[:100])
# print("False alarm prob = {:.20f}".format(ls.false_alarm_probability(power.max())))


a = x[55:90]
b = y[55:90]
c = e[55:90]
print ('phase =',phase[55:90])

# subtract off integer part of MJD
# otherwise floating point errs in minimisation
tfloor = int(a.min())
a -= tfloor

# Fit Gaussian to the min value of dfdt_ts
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
plt.plot(a, gaussian(a, *best_vals), label="fit")
plt.plot(a, g1, label="data")
plt.legend()
plt.show()

mid_ecl = tfloor + best_vals[1]
mid_ecl_err = np.sqrt(covar[1,1])
print("Eclipse Time = ",mid_ecl, "+/-", mid_ecl_err)