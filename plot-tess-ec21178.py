# Program to fit gaussian to light curve- first attempt!
# D Sahman Nov 2018
#=======================================================

import numpy as np
import glob
import matplotlib
from matplotlib import pyplot as plt

# set print option to print to 12 decimal places
# ==============================================
np.set_printoptions(precision=12)

# The lines below are alternative file load methods
# =================================================

# name, mjd, flag, expose, ccd, fwhm, beta, naper, x, y, xm, ym, exm, eym, counts, sigma, sky, nsky, nrej, worst, errorflag, naper2, x2, y2, xm2, ym2, exm2, eym2, counts2, sigma2, sky2, nsky2, nrej2, worst2, errorflag2 = np.loadtxt('oycar_1.log', unpack=True)
# data = np.loadtxt('oycar_1.log', usecols=[0,1,4,7,14,15,20,27,28,34])

# Load file
# =========
print ('Please enter filename: ')
filename=input()
data = np.loadtxt(filename)

a = data.shape
rows = a[0]
cols = a[1]


# Plot all three CCDs
# ===================

fig, axarr = plt.subplots()
mjd1 = data[:,:1]
counts1 = data[:,1:2]
err1 = data[:,2:3]

axarr.plot(mjd1,counts1, linestyle='-', color='red')
axarr.set_ylabel('counts')
axarr.set_xlabel('MJD')
plt.show()
