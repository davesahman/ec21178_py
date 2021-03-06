#!/usr/bin/env python3

#  Program to plot O-C curve for EC21178
#  Uses output file "eclipse_times.txt" from 
#  program eclipse_measure.py



import numpy as np
from numpy import exp, linspace,random, math
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
import glob
from hipercam.hlog import Hlog
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import gaussian_fwhm_to_sigma

lcfile = "eclipse_times.txt"

x,y,e = np.loadtxt(lcfile, unpack=True, usecols=(0,1,2))

# print(x)
# print(y) 
# print(e)

# plt.plot(x,y)
# plt.scatter(x,y, s=3, marker='o')
# plt.errorbar(x, y, yerr=e, fmt='o')
# plt.show()

fit, covars= np.polyfit(x,y,1, cov=True)
print('fit',fit)
print('covars',covars)
print('covars[0,0]',covars[0,0])
print('covars[1,1]',covars[1,1])
print('np.sqrt(abs(covars[0,0]))',np.sqrt(abs(covars[0,0])))
print('np.sqrt(abs(covars[1,1]))',np.sqrt(abs(covars[1,1])))

y2 = x*fit[0] + fit[1]

'''
plt.scatter(x,y, s=3, marker='o', label='data')
plt.plot(x,y2,label='fit')
plt.legend()
plt.show()
'''
plt.figure(figsize=(15,7))
plt.scatter(x,(y-y2)*86400)
plt.title('EC21178   O-C Plot')
plt.axhline(y=0,color='black')
plt.xlabel('Eclipse No.')
plt.ylabel('Time (sec)')
plt.show()


sigt = np.sqrt(abs(covars[0,0]))
sigp = np.sqrt(covars[1,1])
N = round((-1*sigt*sigp)/(sigp**2))
print('N=',N)
