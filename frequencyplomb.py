# code for Lomb Scargle Periodogram of the lightcurve data
# not essential to run, but the PSD accuracy is higher for the python version.

import matplotlib.pyplot as plt
import numpy as np
import csv
import scipy.signal as signal
import glob
from scipy.io import savemat

from numpy import genfromtxt
filepath='lightcurves/'
filenames = sorted(glob.glob(filepath+'lcvold*.dat')) # loading all lightcurve data filenames

c=0
a= np.zeros((300,1))


for nam in filenames:
    data = genfromtxt(nam, delimiter=',')
    data[:,0]=data[:,0]*60;
    data[:,1]=data[:,1]-np.mean(data[:,1])
    f=np.linspace(7e-5*np.pi, 1.1e-4*np.pi, 1000000)    
    # f=np.linspace(1e-6*np.pi, 1e-4*np.pi, 1000000)    # other options for sampling frequencies
    pg= signal.lombscargle(data[:,0],data[:,1], f)
    f=f/(2*np.pi)
    a[c]=f[np.argmax(pg)]
    c=c+1

# saving the PSD data
savemat("j2fr.mat",{"fj2": a})