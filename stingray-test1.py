#!/usr/bin/env python3
import numpy as np
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
# %matplotlib inline
font_prop = font_manager.FontProperties(size=16)


# Generate light curve object
dt = 0.03125  # seconds
exposure = 8.  # seconds
times = np.arange(0, exposure, dt)  # seconds

signal = 300 * np.sin(2.*np.pi*times/0.5) + 1000  # counts/s
noisy = np.random.poisson(signal*dt)  # counts
lc = Lightcurve(times, noisy)
'''
fig, ax = plt.subplots(1,1,figsize=(10,6))
ax.plot(lc.time, lc.counts, lw=2, color='blue')
ax.set_xlabel("Time (s)", fontproperties=font_prop)
ax.set_ylabel("Counts (cts)", fontproperties=font_prop)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.tick_params(which='major', width=1.5, length=7)
ax.tick_params(which='minor', width=1.5, length=4)
plt.show()

ps = Powerspectrum(lc)
print(ps)
print("\nSize of positive Fourier frequencies:", len(ps.freq))
print("Number of data points per segment:", ps.n)

print(ps.freq) # Array of mid-bin frequencies that the Fourier transform samples.
print(ps.power) # Power spectrum
print(ps.df) # The frequency resolution.
print(ps.m) # Number of power spectra averaged together
print(ps.n) # The number of data points (time bins) in one segment of the light curve
print(ps.nphots1) # The total number of photons in the light curve

fig, ax1 = plt.subplots(1,1,figsize=(9,6), sharex=True)
ax1.plot(ps.freq, ps.power, lw=2, color='blue')
ax1.set_ylabel("Frequency (Hz)", fontproperties=font_prop)
ax1.set_ylabel("Power (raw)", fontproperties=font_prop)
ax1.set_yscale('log')
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.tick_params(which='major', width=1.5, length=7)
ax1.tick_params(which='minor', width=1.5, length=4)
for axis in ['top', 'bottom', 'left', 'right']:
    ax1.spines[axis].set_linewidth(1.5)
plt.show()
'''

long_dt = 0.03125  # seconds
long_exposure = 1600.  # seconds
long_times = np.arange(0, long_exposure, long_dt)  # seconds

# In count rate units here
long_signal = 300 * np.sin(2.*np.pi*long_times/0.5) + 1000

# Multiply by dt to get count units, then add Poisson noise
long_noisy = np.random.poisson(long_signal*dt)

long_lc = Lightcurve(long_times, long_noisy)
'''
fig, ax = plt.subplots(1,1,figsize=(10,6))
ax.plot(long_lc.time, long_lc.counts, lw=2, color='blue')
ax.set_xlim(0,20)
ax.set_xlabel("Time (s)", fontproperties=font_prop)
ax.set_ylabel("Counts (cts)", fontproperties=font_prop)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.tick_params(which='major', width=1.5, length=7)
ax.tick_params(which='minor', width=1.5, length=4)
plt.show()
'''

avg_ps = AveragedPowerspectrum(long_lc, 8.)

fig, ax1 = plt.subplots(1,1,figsize=(9,6))
ax1.plot(avg_ps.freq, avg_ps.power, lw=2, color='blue')
ax1.set_xlabel("Frequency (Hz)", fontproperties=font_prop)
ax1.set_ylabel("Power (raw)", fontproperties=font_prop)
ax1.set_yscale('log')
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.tick_params(which='major', width=1.5, length=7)
ax1.tick_params(which='minor', width=1.5, length=4)
for axis in ['top', 'bottom', 'left', 'right']:
    ax1.spines[axis].set_linewidth(1.5)
plt.show()