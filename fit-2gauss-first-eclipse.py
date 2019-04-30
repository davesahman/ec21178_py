#!/usr/bin/env python3

import numpy as np
from numpy import exp, linspace,random
import scipy
from scipy.optimize import curve_fit
import glob
import matplotlib
from matplotlib import pyplot as plt
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

MJD        = 58320
Title      = 'EC21178'
Norm_start = 0.701
device     = '1/xs'
period     = True
lcfile="../data/EC21178-5417.txt"

x,y,e = np.loadtxt(lcfile, unpack=True, usecols=(0,1,2))
'''
data = np.loadtxt(filename)
a = data.shape
rows = a[0]
cols = a[1]

mjd1 = data[:,:1]
counts1 = data[:,1:2]
err1 = data[:,2:3]
'''

x = x[:140]
y = y[:140]
e = e[:140]

# Calculate derivative of the signal
# ==================================
smooth_counts = convolve(y, Box1DKernel(10))
df = smooth_counts[1:] - smooth_counts[:-1]
dt = x[1:] - x[:-1]
dfdt = df/dt
mid_mjd = 0.5 * x[:-1] + 0.5 * x[1:]

# subtract off integer part of MJD
# otherwise floating point errs in minimisation
tfloor = int(mid_mjd.min())
print("tfloor =",tfloor)
mid_mjd -= tfloor
x -= tfloor
# trim the data to exclude noisy end points
mid_mjd = mid_mjd[5:-6]
dfdt = dfdt[5:-6]
# Smooth the derivative
dfdt = convolve(dfdt, Box1DKernel(10))

# Fit two Gaussians to the max and min values of dfdt_ts
#-------------------------------------------------------
# find min and max values of dfdt_ts
min_arg = np.argmin(dfdt)
t_min = mid_mjd[min_arg]
max_arg = np.argmax(dfdt)
t_max = mid_mjd[max_arg]


# get two gaussians around peak to fit
side = 5
t1 = mid_mjd[min_arg-side:min_arg+side]
t2 = mid_mjd[max_arg-side:max_arg+side]
g1 = np.fabs(dfdt[min_arg-side:min_arg+side])
g2 = np.fabs(dfdt[max_arg-side:max_arg+side])

# First gaussian fit
# Initial guesses
amp = -1 * dfdt[min_arg]
cen = t_min
wid = est_fwhm(t1, g1) * gaussian_fwhm_to_sigma
init_vals = [amp ,cen, wid]
best_vals, covar = curve_fit(gaussian, t1, g1, p0=init_vals)
print("Start of eclipse")
print("================")
print("init_vals = ", init_vals)
print ("best vals = ", best_vals)
print(" amp1 = %.3f +/- %.3f" % (best_vals[0], np.sqrt(covar[0,0])))
print(" cen1 = %.3f +/- %.3f" % (best_vals[1], np.sqrt(covar[1,1])))
print(" wid1 = %.3f +/- %.3f" % (best_vals[2], np.sqrt(covar[2,2])))

# Second Gaussain fit

# Initial guesses
amp2 = dfdt[max_arg]
cen2 = t_max
wid2 = est_fwhm(t2, g2) * gaussian_fwhm_to_sigma
print("cen2 = ", cen2)
init_vals2 = [amp2 ,cen2, wid2]

best_vals2, covar2 = curve_fit(gaussian, t2, g2, p0=init_vals2)
print("End of eclipse")
print("==============")
print("init_vals2 = ", init_vals2)
print ("best vals2 = ", best_vals2)
print(" amp2 = %.3f +/- %.3f" % (best_vals2[0], np.sqrt(covar2[0,0])))
print(" cen2 = %.3f +/- %.3f" % (best_vals2[1], np.sqrt(covar2[1,1])))
print(" wid2 = %.3f +/- %.3f" % (best_vals2[2], np.sqrt(covar2[2,2])))

fig, axs = plt.subplots(2,2,figsize=(10,7))
axs[0,0].plot(x, y, linestyle='-', color='red',marker='o')
ax = axs[0, 0]
# ax.set_xlabel("Time (days)")
ax.set_title(" Light curve - EC21178")
axs[0, 1].plot(mid_mjd, dfdt, linestyle='-', color ='red',marker='o')
ax = axs[0, 1]
# ax.set_xlabel("Time (days)")
ax.set_title("Smoothed Derivative")
axs[1, 0].plot(t1, gaussian(t1, *best_vals), label="fit")
axs[1, 0].plot(t1, g1, label="data")
axs[1, 0].legend()
ax = axs[1, 0]
ax.set_xlabel("Time (MJD - {})".format(tfloor))
ax.set_title("First Gaussian Fit")
axs[1, 1].plot(t2, gaussian(t2, *best_vals2), label="fit")
axs[1, 1].plot(t2, g2, label="data")
axs[1, 1].legend()
ax = axs[1, 1]
ax.set_xlabel("Time (MJD - {})".format(tfloor))
ax.set_title("Second Gaussian Fit")
plt.show()





'''
# Plot all three CCDs
# ===================
fig, axarr = plt.subplots()
axarr.plot(x,y, linestyle='-', color='red')
axarr.set_ylabel('counts')
axarr.set_xlabel('MJD')
plt.show()
'''
exit()