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



