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
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum
from astropy.timeseries import LombScargle
from astropy.time import Time
from scipy.optimize import curve_fit
from scipy import interpolate
from hipercam.hlog import Hlog, Tseries
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from scipy.stats import pearsonr

np.set_printoptions(precision=12)

def gaussian(x, amp, cen, wid):
    """
    Set up Gaussian model
    """
    return amp * exp (-(x-cen)**2/(2*wid**2))

def est_fwhm(x, y):
    """
    Estimate FWHM of a Gaussian
    """
    half_max = 0.5*y.max()
    within_halfmax = y > half_max
    x_within_halfmax = x[within_halfmax]
    return x_within_halfmax.max() - x_within_halfmax.min()

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
# x *= 1440. 

# Plot raw light curve

plt.figure(figsize=(16,8))
plt.plot(x,y)
plt.title('EC21178 TESS Raw Light Curve')
plt.xlabel('Time (BJD-TDB)')
plt.ylabel('FLux')
plt.show()


# Find the Average Power Spectrum using Stingray
'''
lc = Lightcurve(x,y)
ps = AveragedPowerspectrum(lc,1.)
#lin_rb_ps = ps.rebin(0.01, method='mean')
fig, ax1 = plt.subplots(1,1,figsize=(9,6), sharex=True)
ax1.plot(ps.freq, ps.power, lw=2, color='blue')
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("Power (raw)")
# ax1.set_yscale('log')
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.tick_params(which='major', width=1.5, length=7)
ax1.tick_params(which='minor', width=1.5, length=4)
for axis in ['top', 'bottom', 'left', 'right']:
    ax1.spines[axis].set_linewidth(1.5)
plt.show()

print(ps.freq)
print(ps.power)
print(ps.df)
print(ps.m)
print(ps.n)
print(ps.nphots1)
'''

# Find an approximate period using a
# Lomb-Scargle periodogram

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
    plt.figure(figsize=(16,8))
    plt.plot(freq, power)
    plt.show()
    '''

period_time = 1/fmax # derive the period
phase = fmax*x # calculate the phase of each data point
phase = np.mod(phase,1) # convert the phase to lie between 0 and 1

tfloor = period_time * int(x[0]/period_time) 
x -= tfloor # subtract integer of first period 
cycle1 = np.floor_divide(x,period_time) # calculate which cycle each data point lies in
cycle_vals = np.unique(cycle1) #  calculate the number of unique cycle values ie number of  eclipses

# loop over each cycle and fit gaussian

ecl = np.zeros([int(cycle_vals.max()),5],dtype=float) # set up ecl array and fill with zeroes
cumy = 0
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
      print("y within 0.4",len(y_fwtm))
      cumy +=len(y_fwtm)
      avgy = cumy/i
    # Gaussian fit
    # Initial guesses
      amp = yr[min_arg]
      cen = t_min
      wid = est_fwhm(x_fwtm, y_fwtm) * gaussian_fwhm_to_sigma
      init_vals = [amp ,cen, wid]
      best_vals, covar = curve_fit(gaussian, x_fwtm, y_fwtm, p0=init_vals)

    # Plot results of gaussian fits
      '''
      plt.figure(figsize=(16,8))
      plt.plot(x_range, gaussian(x_range, *best_vals), label="fit")
      plt.plot(x_range, yr, label="data")
      plt.legend()
      plt.show()
      
      if i > 10:
        print('average number of y values used =', avgy)
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
s = ecl.shape[0]

#Save the ecl array
np.savetxt("eclipse_times.txt",ecl)

# fit line to eclipse times

z = np.polyfit(ecl[:,0], ecl[:,1], 1, cov=True)
p = np.poly1d(z[0])

print('z[0]=',z[0])
print('z[1]=',z[1])
print('z=',z)


# Calculate zeropoint eclipse number using Stu's method

c0 = np.rint(-z[1][0,1]/z[1][0,0])
ephem = p(c0) + tfloor
t0 = p(1) + tfloor
epheme = np.sqrt(z[1][1,1])
periode = np.sqrt(z[1][0,0])

# Print key results

print('Optimal eclipse number          = ',c0)
print('Ephemeris for publication       = ',ephem)
print('Error Ephemeris for publn       = ',epheme)
print('Derived Time of first eclipse   = ',t0)
print('Time of ecl 0                   = ',ecl[0,1]+tfloor)
print('Shape of ecl array              = ',s)
print('Period from linear fit (days)   = ',z[0][0])
print('Error on Period from fit(days)   = ',periode)
print('Lomb Scargle period (days)      = ',period_time)

# Plot O-C curve

plt.figure(figsize=(16,8))
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
ts2 = ts.fold(z[0][0],t0)
ts3 = ts2.bin(1000,'mean')
print('ts3.ye',ts3.ye)
'''
plt.figure(figsize=(10,6))
plt.plot(ts3.t,ts3.ye)
# plt.plot(ts3.t,ts3.y)
#plt.plot(ts3.t-0.5,ts3.y)
#plt.plot(ts3.t+0.5,ts3.y)
# plt.xlim(-0.5,1.5)
plt.title('EC21178 Y-error Binned Light Curve')
plt.xlabel('Phase')
plt.ylabel('FLux')
plt.show()
'''

f = interpolate.interp1d(ts3.t,ts3.y,kind='linear',fill_value="extrapolate")
tsnew = ts.t.copy()
ph = np.mod(((tsnew + (z[0][0]/2) - t0)) / z[0][0],1) - 0.5
ts.y -= f(ph)


ls    = LombScargle(ts.t,ts.y)                   # Create periodogram
fmax  = 30.                                      # Set upper frequency (cycles/min) limit 
nfreq = int(1000*fmax*(ts.t.max()-ts.t.min()))   # Calculate number of frequency steps, oversample x10
freq  = np.linspace(fmax/nfreq,fmax,nfreq)       # Create frequency array
#freq = np.linspace(0.1,5,100)
power = ls.power(freq)                           # Calculate periodogram powers            
fmax  = freq[power==power.max()]                 # Calculate peak in periodogra

# Plot Periodogram

plt.figure(figsize=(16,8))
plt.plot(freq, power)
plt.show()


# Plot noise curve

corr, _ = pearsonr(ecl[:,0], ecl[:,3])
plt.figure(figsize=(16,8))
plt.plot(ts.t,ts.y)
plt.text(1, 80, 'Pearsons Coefficient= (%a)'%(corr))
plt.title('Residual Noise')
plt.xlabel('Time')
plt.ylabel('Counts')
plt.show()


# Subtract linear fit from noise

z = np.polyfit(ts.t, ts.y, 1)
p = np.poly1d(z)
ts.y -= p(ts.t)
ts2 = ts.fold(z[0],t0) 
ts3 = ts2.bin(40,'mean') 

# Plot the folded noise light curve - repeat two periods for clarity

plt.figure(figsize=(16,8))
plt.plot(ts3.t,ts3.y)
plt.plot(ts3.t+1.0,ts3.y)
plt.xlim(-0.5,1.5)
plt.title('Binned and Folded Residual Noise')
plt.xlabel('Phase')
plt.ylabel('Counts')
plt.show()

# Plot eclipse Amplitude and Depths to see if any trend

corr1, _ = pearsonr(ecl[:,0], ecl[:,3])
corr2, _ = pearsonr(ecl[:,0], ecl[:,4])

fig = plt.figure(figsize=(16,8))
ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4])
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.4])

ax1.plot(ecl[:,0],ecl[:,3], '-', marker='o')
ax1.text(.5,.9,'Amplitude (FWHM)',
        horizontalalignment='center',
        transform=ax1.transAxes)
ax1.text(.8,.85,'Pearsons coefficient= (%a)'%(corr1),
        horizontalalignment='center',
        transform=ax1.transAxes)       
ax2.plot(ecl[:,0],ecl[:,4], '-', marker='o')
ax2.text(.5,.9,'Depth',
        horizontalalignment='center',
        transform=ax2.transAxes)
ax2.text(.8,.85,'Pearsons coefficient= (%a)'%(corr2),
        horizontalalignment='center',
        transform=ax2.transAxes)          
ax2.set_xlabel('eclipse number')        
plt.show()

# Plot Amplitude vs Depth and fit linear regression 
# Show Pearsons correlation coefficient

q = np.polyfit(ecl[:,3], ecl[:,4], 1)
qt0 = str(q[0])
p = np.poly1d(q)
fit = p(ecl[:,3])
corr3, _ = pearsonr(ecl[:,3], ecl[:,4])

plt.figure(figsize=(16,8))
plt.scatter(ecl[:,3], ecl[:,4], marker='o')
plt.plot(ecl[:,3], fit, color='r', label='fit')
plt.title('Scatter Plot of Amplitude versus Depth')
plt.text(300, 0.014, 'Slope of Fit = (%a)'%(qt0))
plt.text(300, 0.0135, 'Pearsons Coefficient= (%a)'%(corr3))
plt.xlabel('Amplitude')
plt.ylabel('FWHM')
plt.show()

