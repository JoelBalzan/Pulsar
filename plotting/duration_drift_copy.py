import glob
import os
import sys

import lmfit
import matplotlib.pyplot as plt
import numpy as np
import psrchive
from lmfit import Model
from lmfit.lineshapes import lorentzian
from lmfit.model import save_modelresult
from matplotlib import gridspec
from matplotlib.patches import Ellipse
from scipy import signal
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_widths

duration_intra_drift = np.load('duration_intra_drift.npy')
duration_intra_drift = np.array(duration_intra_drift)
duration_sub_drift = np.load('duration_sub_drift.npy')
duration_sub_drift = np.array(duration_sub_drift)


duration_intra_drift_pos = []
duration_intra_drift_neg = []
for i in range(np.shape(duration_intra_drift)[0]):
    if duration_intra_drift[i,1] > 0:
        duration_intra_drift_pos.append(duration_intra_drift[i,:])
    else:
        duration_intra_drift_neg.append(duration_intra_drift[i,:])
duration_intra_drift_pos = np.array(duration_intra_drift_pos)
duration_intra_drift_neg = np.array(duration_intra_drift_neg)

duration_sub_drift_pos = []
duration_sub_drift_neg = []
for i in range(np.shape(duration_sub_drift)[0]):
    if duration_sub_drift[i,1] > 0:
        duration_sub_drift_pos.append(duration_sub_drift[i,:])
    else:
        duration_sub_drift_neg.append(duration_sub_drift[i,:])
duration_sub_drift_pos = np.array(duration_sub_drift_pos)
duration_sub_drift_neg = np.array(duration_sub_drift_neg)

duration_intra_drift_neg_new = []
for i in range(np.shape(duration_intra_drift_neg)[0]):
    if -6 < duration_intra_drift_neg[i,1] < 0:
        duration_intra_drift_neg_new.append(duration_intra_drift_neg[i,:])
duration_intra_drift_neg_new = np.array(duration_intra_drift_neg_new)
### PLOTTING ###
A4x, A4y = 8.27, 11.69

fig = plt.figure(figsize=(A4y/2, A4x/2), dpi=600)
plt.plot(duration_intra_drift_neg_new[:,0], duration_intra_drift_neg_new[:,1], 'o', markersize=2, color='green')
#plt.plot(duration_intra_drift_pos[:,0], duration_intra_drift_pos[:,1], 'o', markersize=2)
slope1, intercept1 = np.round(np.polyfit(duration_intra_drift_neg_new[:,0], duration_intra_drift_neg_new[:,1], 1), 2)
mb = np.round(np.polyfit(duration_intra_drift_neg_new[:,0], duration_intra_drift_neg_new[:,1], 2), 2)
plt.plot(np.unique(duration_intra_drift_neg_new[:,0]), np.poly1d(np.polyfit(duration_intra_drift_neg_new[:,0], duration_intra_drift_neg_new[:,1], 2))(np.unique(duration_intra_drift_neg_new[:,0])), linestyle='--',
         label = r'$\frac{dt}{d\nu}$ = %sDur$^{2}_{tot}$ + %sDur$_{tot}$ - %s'%(mb[0], mb[1], np.abs(mb[2])), 
         color='orange'
         )
plt.plot(np.unique(duration_intra_drift_neg_new[:,0]), np.poly1d(np.polyfit(duration_intra_drift_neg_new[:,0], duration_intra_drift_neg_new[:,1], 1))(np.unique(duration_intra_drift_neg_new[:,0])), linestyle='-',
         label = r'$\frac{dt}{d\nu}$ = %sDur$_{tot}$ - %s'%(slope1, np.abs(intercept1)), 
         color='red'
         )
#plt.plot(np.unique(duration_intra_drift_pos[:,0]), np.poly1d(np.polyfit(duration_intra_drift_pos[:,0], duration_intra_drift_pos[:,1], 1))(np.unique(duration_intra_drift_pos[:,0])), linestyle='-')
#plt.xscale('log')
#plt.yscale('log')
plt.ylim(-5, -4.2)
#plt.xlim(0,2.5)
plt.xlabel('Duration (ms)')
plt.ylabel('Drift Rate (ms/MHz)')
plt.legend()
plt.savefig('duration_intra_drift.png', dpi=600, bbox_inches='tight')
print('duration_intra_drift.png')


fig = plt.figure(figsize=(A4y/2, A4x/2), dpi=600)
plt.plot(duration_sub_drift_neg[:,0], duration_sub_drift_neg[:,1], 'o', markersize=2)
plt.plot(duration_sub_drift_pos[:,0], duration_sub_drift_pos[:,1], 'o', markersize=2)
slope, intercept = np.round(np.polyfit(duration_sub_drift_neg[:,0], duration_sub_drift_neg[:,1], 1), 2)
plt.plot(np.unique(duration_sub_drift_neg[:,0]), np.poly1d(np.polyfit(duration_sub_drift_neg[:,0], duration_sub_drift_neg[:,1], 1))(np.unique(duration_sub_drift_neg[:,0])), linestyle='-',
         #label = r'$\frac{d\nu}{dt}_neg$ = %sDur$_{tot}$ - %s'%(slope, np.abs(intercept)), 
         color='red',
         label = r'$\frac{d\nu}{dt}_{sad}$')
slope, intercept = np.round(np.polyfit(duration_sub_drift_pos[:,0], duration_sub_drift_pos[:,1], 1), 2)
plt.plot(np.unique(duration_sub_drift_pos[:,0]), np.poly1d(np.polyfit(duration_sub_drift_pos[:,0], duration_sub_drift_pos[:,1], 1))(np.unique(duration_sub_drift_pos[:,0])), linestyle='-',
         #label = r'$\frac{d\nu}{dt}_pos$ = %sDur$_{tot}$ + %s'%(slope, intercept), 
         color='green',
         label = r'$\frac{d\nu}{dt}_{happy}$')

#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Duration (ms)', fontsize=10)
plt.ylabel('Sub-Burst Drift Rate (MHz/ms)', fontsize=10)
plt.legend()
plt.savefig('duration_sub_drift.pdf', dpi=600, bbox_inches='tight')
print('duration_sub_drift.png')