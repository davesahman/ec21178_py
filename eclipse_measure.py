#!/usr/bin/env python3

# This program analyses the TESS light curve of EC 21178-5417
# Author D Sahman October 2019

import sys
import numpy as np
import astropy
import scipy
import glob
import hipercam as hcam
import matplotlib.pyplot as plt
from numpy import exp, linspace, random, math
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit
from hipercam.hlog import Hlog
from hipercam.hlog import Tseries
from scipy import interpolate
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
    """
    Finds the element in the array `a` closest to the scalar value `a0`
    """
    idx = np.abs(a - a0).argmin()
    return a.flat[idx]

MJD        = 58324
lcfile     = '../data/EC21178-5417.txt'
Title      = 'EC21178'
Norm_start = 0.701
device     = '1/xs'
period     = True

x,y,e = np.loadtxt(lcfile, unpack=True, usecols=(0,1,2))

if period:
    # x *= 1440.                                # Change to minutes
    ls    = LombScargle(x,y,e)                  # Create periodogram
    fmax  = 30.                                 # Set upper frequency (cycles/min) limit 
    nfreq = int(1000*fmax*(x.max()-x.min()))    # Calculate number of frequency steps, oversample x10
    freq  = np.linspace(fmax/nfreq,fmax,nfreq)  # Create frequency array
    #freq = np.linspace(0.1,5,100)
    power = ls.power(freq)                      # Calculate periodogram powers            
    fmax  = freq[power==power.max()]            # Calculate peak in periodogram

    # Plot Periodogram
    '''
    plt.figure(figsize=(20,10))
    plt.plot(freq, power)
    plt.show()
    '''
period_time = 1/fmax
phase = fmax*x
phase = np.mod(phase,1)

tfloor = period_time * int(x[0]/period_time)
x -= tfloor # subtract integer of first period 
cycle1 = np.floor_divide(x,period_time)
cycle_vals = np.unique(cycle1)

# loop over each cycle and fit gaussian

ecl = np.zeros([int(cycle_vals.max()),5],dtype=float) # set up ecl array and fill with zeroes

i = 0
while i < (cycle_vals.max()-1):
    i += 1
    if i in cycle1:
      if i == 33 or i == 86 or i ==93 or i == 146 or i ==149 or i==152 or i ==153 or i == 155:  # these eclipses are bad
         continue
      c_range = cycle1 == i
      x_range = x[c_range]
      y_range = y[c_range]
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
      plt.figure(figsize=(20,10))
      plt.plot(x_fwtm, gaussian(x_fwtm, *best_vals), label="fit")
      plt.plot(x_range, yr, label="data")
      plt.legend()
      plt.show()
      '''
      
      ecl[i,0] +=i # number of eclipse
      ecl[i,1] +=best_vals[1] # centre of eclipse - need to add tfloor for true time
      ecl[i,2] +=np.sqrt(covar[1,1]) # error on centre
      ecl[i,3] +=best_vals[0] # amplitude
      ecl[i,4] += best_vals[2] # width

    # scale y values to y.min
      '''
      print('y[c_range]',y[c_range])
      print('y.min()',y.min())
      print('y[c_range].min()',y[c_range].min())
      
      y[c_range] -= (y[c_range].min() - y.min()) 
      '''

# remove rows in ecl array with all zero entries
ecl = ecl[~np.all(ecl==0, axis=1)]

#Save the ecl array
np.savetxt("eclipse_times.txt",ecl)

# fit line to eclipse times

z = np.polyfit(ecl[:,0], ecl[:,1], 1)
p = np.poly1d(z)
t0 = z[1] + z[0] + tfloor # time of first eclipse using fit
s = ecl.shape[0]

# Print results of fit

print('z',z)
print('shape of ecl array = ',s)
print('period from linear fit (days)          = ',z[0])
print('Lomb Scargle period (days) = ',period_time)
print('time of ecl 0',ecl[0,1])
print('tfloor = ',tfloor)
print('Time of first eclipse',t0)

# Plot O-C curve

plt.figure(figsize=(20,10))
plt.scatter(ecl[:,0],(ecl[:,1]-p(ecl[:,0]))*86400)
plt.axhline(y=0,color='black')
plt.title('EC21178   O-C Plot')
plt.xlabel('Eclipse No.')
plt.ylabel('Time (sec)')
plt.show()

# Measure Noise by subtracting binned average light curve 
# Create Hipercam Tseries object to do folding

mask = np.zeros_like(x).astype('int')
ts = Tseries(x,y,e,mask)
ts2 = ts.fold(z[0],t0)
ts3 = ts2.bin(20,'mean')
f = interpolate.interp1d(ts3.t,ts3.y,kind='linear',fill_value="extrapolate")
tsnew = ts.t.copy()
ph = np.mod(((tsnew + (z[0]/2) - t0)) / z[0],1) - 0.5
ts.y -= f(ph)

# Plot noise curve - repeat two periods for clarity
'''
plt.figure(figsize=(16,8))
plt.plot(ts.t,ts.y)
plt.title('Residual Noise')plt.title('EC21178   O-C Plot')
plt.xlabel('Time')
plt.ylabel('Counts')

plt.show()
'''
# Subtract linear fit from noise

z = np.polyfit(ts.t, ts.y, 1)
p = np.poly1d(z)
ts.y -= p(ts.t)
ts2 = ts.fold(z[0],t0) 
ts3 = ts2.bin(40,'mean') 

# Plot the folded noise light curve
'''
plt.figure(figsize=(12,6))
plt.plot(ts3.t,ts3.y)
plt.plot(ts3.t+1.0,ts3.y)
plt.xlim(-0.5,1.5)
plt.title('Binned and Folded Residual Noise')
plt.xlabel('Phase')
plt.ylabel('Counts')

plt.show()
'''
# Plot eclipse Amplitude and Depths to see if any trend
'''
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4])
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.4])

ax1.plot(ecl[:,0],ecl[:,3], '-', marker='o')
ax1.text(.5,.9,'Amplitude (FWHM)',
        horizontalalignment='center',
        transform=ax1.transAxes)
ax2.plot(ecl[:,0],ecl[:,4], '-', marker='o')
ax2.text(.5,.9,'Depth',
        horizontalalignment='center',
        transform=ax2.transAxes)
ax2.set_xlabel('eclipse number')        
plt.show()

# Plot Amplitude vs Depth and fit linear regression to see if 
# there is a correlation 

q = np.polyfit(ecl[:,3], ecl[:,4], 1)
qt0 = str(q[0])
p = np.poly1d(q)
fit = p(ecl[:,3])

plt.figure(figsize=(20,10))
plt.scatter(ecl[:,3], ecl[:,4], marker='o')
plt.plot(ecl[:,3], fit, color='r', label='fit')
plt.title('Scatter Plot of Amplitude versus Depth')
plt.text(300, 0.014, 'Slope of Fit = (%a)'%(qt0))
plt.xlabel('Amplitude')
plt.ylabel('FWHM')
plt.show()
'''