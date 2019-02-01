#!/usr/bin/env python

import numpy as np
from ppgplot import *
from astropy.stats import LombScargle

MJD        = 58466
lcfile     = 'run014.dat'
Title      = 'WDJ044832.11-105349.85'
Norm_start = 0.701
device     = '1/xs'
period     = True

x,y,e = np.loadtxt(lcfile, unpack=True, usecols=(0,2,3))

x -= MJD

if len(x[x<Norm_start])>0:
    e /= y[x<Norm_start].mean()
    y /= y[x<Norm_start].mean()

pgopen(device)

if device != '1/xs':
    pgscr(0,1,1,1)
    pgscr(1,0,0,0)

pgscf(2)
pgsch(1.6)
pgslw(4)

pgenv(x.min(),x.max(),y.min()*0.99,y.max()*1.01)
#pgenv(x.min(),x.max(),0.98,1.02)
pglab('MJD - '+str(MJD),'Normalised counts',Title)

pgslw(3)
pgsch(1.2)
pgpt(x,y,17)
pgerry(x,y+e,y-e,0)

pgclos()

if period:
    x *= 1440.                                 # Change to minutes
    ls    = LombScargle(x,y,e)                 # Create periodogram
    fmax  = 1.                                 # Set upper frequency (cycles/min) limit 
    nfreq = int(10*fmax*(x.max()-x.min()))     # Calculate number of frequency steps, oversample x10
    freq  = np.linspace(fmax/nfreq,fmax,nfreq) # Create frequency array
    power = ls.power(freq)                     # Calculate periodogram powers
    fmax  = freq[power==power.max()][0]        # Calculate peak in periodogram
    print("Peak in periodogram at {:5.3f} cycles / min. Period of {:5.3f} minutes".format(fmax,1./fmax))

    pgopen('2/xs')
    pgscf(2)
    pgsch(1.6)
    pgslw(4)

    pgenv(freq.min(),freq.max(),0,power.max())
    pglab('Cycles/min','Power','')

    pgline(freq,power)
    pgclos()

