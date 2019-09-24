#!/usr/bin/env python3

# =========================================================
# Program to calculate average light curve and subtract to
# show underlying trends. Uses file 'ec21178.log' generated
# by eclipse_measure.py
# =========================================================

import sys
import numpy as np
import matplotlib.pyplot as plt
import hipercam as hcam
from hipercam.hlog import Tseries
from scipy import interpolate

x,y,e = np.loadtxt('ec21178.log', unpack=True, usecols=(0,1,2))
mask = np.zeros_like(x).astype('int')
ts = Tseries(x,y,e,mask)
period = 0.15452941957207128
t1 = 0.196946800588
ts2 = ts.fold(period,t1)
ts3 = ts2.bin(50,'mean')
f = interpolate.interp1d(ts3.t,ts3.y,kind='linear',fill_value="extrapolate")
tsnew = ts.t.copy()
ph = np.mod(((tsnew + (period/2) - t1)) / period,1) - 0.5
# sys.exit() 
ynew = f(ph) 
ydiff = ts.y - ynew
print('tsnew',tsnew)
print('ph=',ph)
print('ynew=',ynew)
sys.exit()
plt.figure(figsize=(16,8))
plt.plot(ts.t,ts.y)
# plt.plot(ts.t,ynew)
# plt.plot(ts.t,ydiff)

plt.show()

dat = np.stack((ts.t,ydiff,e),axis=-1)
np.savetxt("ydiff.log",dat)