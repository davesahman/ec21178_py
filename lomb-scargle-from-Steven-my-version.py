#!/usr/bin/env python3

import numpy as np
from astropy.stats import LombScargle
import matplotlib.pyplot as plt

MJD        = 58324
lcfile     = '../data/EC21178-5417.txt'
Title      = 'WDJ044832.11-105349.85'
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
   
print("False alarm prob = {:.20f}".format(ls.false_alarm_probability(power.max())))

fig,axs = plt.subplots(1,2,figsize=(10,7))
axs[0].plot(freq,power)
ax =axs[0]
ax.set_xlabel("Freq ")
ax.set_ylabel("Power")
ax.set_title("Lomb-Scargle periodogram")

axs[1].plot(phase,y,linewidth=0.2)
ax =axs[1]
ax.set_xlabel("Phase")
ax.set_ylabel("Counts")
ax.set_title("Phase folded light curve")

plt.show()