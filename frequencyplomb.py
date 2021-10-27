import matplotlib.pyplot as plt
import numpy as np
import csv
import scipy.signal as signal
import glob
from scipy.io import savemat

from numpy import genfromtxt
filepath='lightcurves/'
filenames = sorted(glob.glob(filepath+'lcvold*.dat'))
# a = np.zeros((300,1000000))

c=0

# data = genfromtxt('lightcurves/lcvold005.dat', delimiter=',')
# data[:,0]=data[:,0]*60;
# data[:,1]=data[:,1]-np.mean(data[:,1])
# f=np.linspace(8e-5*np.pi, 2e-4*np.pi, 1000000)
# pg=signal.lombscargle(data[:,0],data[:,1], f)
# f=f/(2*np.pi)
# a=f[np.argmax(pg)]
# print(a)
a= np.zeros((300,1))


for nam in filenames:
    data = genfromtxt(nam, delimiter=',')
    # data = genfromtxt('lightcurves/lcvold300.dat', delimiter=',')
    data[:,0]=data[:,0]*60;
    data[:,1]=data[:,1]-np.mean(data[:,1])
    f=np.linspace(7e-5*np.pi, 1.1e-4*np.pi, 1000000)
    # f=np.linspace(1e-6*np.pi, 1e-4*np.pi, 1000000)
    pg=signal.lombscargle(data[:,0],data[:,1], f)
    f=f/(2*np.pi)
    # a[c,:]=pg
    a[c]=f[np.argmax(pg)]
    c=c+1
    print(c)
# print(a)
savemat("j2fr.mat",{"fj2": a})
# my_data = genfromtxt('lightcurves/lcvold005.dat', delimiter=',')
# print(my_data[:,0])

# plt.plot(f,pg)
# plt.show()



# peaks,_=signal.find_peaks(pg)
# print(filenames[1])

# plt.plot(f,pg)
# plt.show()

# rng = np.random.default_rng()
# A = 2.
# w = 1.
# phi = 0.5 * np.pi
# nin = 1000
# nout = 100000
# frac_points = 0.9  # Fraction of points to select
# x = np.linspace(0.01, 10*np.pi, nin)
# r = rng.standard_normal(nin)
# x = x[r >= frac_points]

# y = A * np.sin(w*x+phi)
# f = np.linspace(0.01, 10, nout)
# import scipy.signal as signal
# pgram = signal.lombscargle(x, y, f, normalize=True)
# plt.subplot(2, 1, 1)
# plt.plot(x, y, 'b+')
# plt.subplot(2, 1, 2)
# plt.plot(f, pgram)
# plt.show()
