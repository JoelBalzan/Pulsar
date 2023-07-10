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

duration_intra_drift = np.delete(duration_intra_drift, np.where(np.abs(duration_intra_drift[:,3]/duration_intra_drift[:,2]) > 0.2)[0], 0)
duration_intra_drift = np.delete(duration_intra_drift, np.where(np.abs(duration_intra_drift[:,2]) > 0.0003)[0], 0)


duration_intra_drift_pos = []
duration_intra_drift_neg = []
for i in range(np.shape(duration_intra_drift)[0]):
    if duration_intra_drift[i,2] > 0:
        duration_intra_drift_pos.append(duration_intra_drift[i,:])
    else:
        duration_intra_drift_neg.append(duration_intra_drift[i,:])
duration_intra_drift_pos = np.array(duration_intra_drift_pos)
duration_intra_drift_neg = np.array(duration_intra_drift_neg)
print(np.min(duration_intra_drift_neg[:,2]), np.max(duration_intra_drift_neg[:,2]))
print(np.min(duration_intra_drift_pos[:,2]), np.max(duration_intra_drift_pos[:,2]))
print(duration_intra_drift_pos[np.where(duration_intra_drift_pos[:,2] == np.min(duration_intra_drift_pos[:,2]))[0]])

duration_sub_drift_pos = []
duration_sub_drift_neg = []
for i in range(np.shape(duration_sub_drift)[0]):
    if duration_sub_drift[i,2] > 0:
        duration_sub_drift_pos.append(duration_sub_drift[i,:])
    else:
        duration_sub_drift_neg.append(duration_sub_drift[i,:])
duration_sub_drift_pos = np.array(duration_sub_drift_pos)
duration_sub_drift_neg = np.array(duration_sub_drift_neg)
print(np.min(duration_sub_drift_neg[:,2]), np.max(duration_sub_drift_neg[:,2]))
print(np.min(duration_sub_drift_pos[:,2]), np.max(duration_sub_drift_pos[:,2]))

#duration_intra_drift_neg_new = []
#for i in range(np.shape(duration_intra_drift_neg)[0]):
#    if -1 < duration_intra_drift_neg[i,1] < 1:
#        duration_intra_drift_neg_new.append(duration_intra_drift_neg[i,:])
#duration_intra_drift_neg_new = np.array(duration_intra_drift_neg_new)


### PLOTTING ###
A4x, A4y = 8.27, 11.69
fig = plt.figure(figsize=(A4y/2, A4x/2), dpi=600)
plt.errorbar(duration_intra_drift_neg[:,0], duration_intra_drift_neg[:,2], xerr=np.abs(duration_intra_drift_neg[:,1]), yerr=duration_intra_drift_neg[:,3], fmt='none', color='gray')
#plt.errorbar(duration_intra_drift_pos[:,0], duration_intra_drift_pos[:,2], xerr=np.abs(duration_intra_drift_pos[:,1]), yerr=duration_intra_drift_pos[:,3], fmt='none', color='gray')
plt.plot(duration_intra_drift_neg[:,0], duration_intra_drift_neg[:,2], 'o', markersize=2, color='purple')
#plt.plot(duration_intra_drift_pos[:,0], duration_intra_drift_pos[:,2], 'o', markersize=2, color='magenta')
slope1, intercept1 = np.round(np.polyfit(duration_intra_drift_neg[:,0], duration_intra_drift_neg[:,2], 1), 10)
plt.plot(np.unique(duration_intra_drift_neg[:,0]), np.poly1d(np.polyfit(duration_intra_drift_neg[:,0], duration_intra_drift_neg[:,2], 1))(np.unique(duration_intra_drift_neg[:,0])), linestyle='-',
         label = r'$\frac{dt}{d\nu}$ = %sDur$_{sub}$ - %s'%(slope1, np.abs(intercept1)), 
         color='red'
         )
#mb = np.round(np.polyfit(duration_intra_drift_neg[:,0], duration_intra_drift_neg[:,2], 2), 4)
#plt.plot(np.unique(duration_intra_drift_neg[:,0]), np.poly1d(np.polyfit(duration_intra_drift_neg[:,0], duration_intra_drift_neg[:,2], 2))(np.unique(duration_intra_drift_neg[:,0])), linestyle='--',
#         label = r'$\frac{dt}{d\nu}$ = %sDur$^{2}_{sub}$ + %sDur$_{sub}$ - %s'%(mb[0], mb[1], np.abs(mb[2])), 
#         color='orange'
#         )
#slope1, intercept1 = np.round(np.polyfit(duration_intra_drift_pos[:,0], duration_intra_drift_pos[:,2], 1), 10)
#plt.plot(np.unique(duration_intra_drift_pos[:,0]), np.poly1d(np.polyfit(duration_intra_drift_pos[:,0], duration_intra_drift_pos[:,2], 1))(np.unique(duration_intra_drift_pos[:,0])), linestyle='-',
#         label = r'$\frac{dt}{d\nu}$ = %sDur$_{sub}$ + %s'%(slope1, np.abs(intercept1)), 
#         color='green'         
#         )
#plt.xscale('log')
#plt.yscale('log')
#plt.ylim(-0.00025, 0.0025)
#plt.ylim(-1, 1)

#plt.xlim(0,2.5)
plt.xlabel('Duration (ms)')
plt.ylabel('Intra-Burst Drift Rate ' + r'$\frac{dt}{d\nu}$' + ' (ms/MHz)')
#plt.legend()
plt.savefig('duration_intra_drift.pdf', dpi=600, bbox_inches='tight')
print('duration_intra_drift.png')




fig = plt.figure(figsize=(A4y/2, A4x/2), dpi=600)
plt.errorbar(duration_sub_drift_neg[:,0], duration_sub_drift_neg[:,2], xerr=np.abs(duration_sub_drift_neg[:,1]), yerr=duration_sub_drift_neg[:,3], fmt='none', color='gray')
plt.errorbar(duration_sub_drift_pos[:,0], duration_sub_drift_pos[:,2], xerr=np.abs(duration_sub_drift_pos[:,1]), yerr=duration_sub_drift_pos[:,3], fmt='none', color='gray')
plt.plot(duration_sub_drift_neg[:,0], duration_sub_drift_neg[:,2], 'o', markersize=2, color='blue')
plt.plot(duration_sub_drift_pos[:,0], duration_sub_drift_pos[:,2], 'o', markersize=2, color='orange')
slope, intercept = np.round(np.polyfit(duration_sub_drift_neg[:,0], duration_sub_drift_neg[:,2], 1), 2)
plt.plot(np.unique(duration_sub_drift_neg[:,0]), np.poly1d(np.polyfit(duration_sub_drift_neg[:,0], duration_sub_drift_neg[:,2], 1))(np.unique(duration_sub_drift_neg[:,0])), linestyle='-',
         label = r'$\frac{d\nu}{dt}_neg$ = %sDur$_{sub}$ - %s'%(slope, np.abs(intercept)), 
         color='red',
         #label = r'$\frac{d\nu}{dt}_{sad}$'
         )
slope, intercept = np.round(np.polyfit(duration_sub_drift_pos[:,0], duration_sub_drift_pos[:,2], 1), 2)
plt.plot(np.unique(duration_sub_drift_pos[:,0]), np.poly1d(np.polyfit(duration_sub_drift_pos[:,0], duration_sub_drift_pos[:,2], 1))(np.unique(duration_sub_drift_pos[:,0])), linestyle='-',
         label = r'$\frac{d\nu}{dt}_pos$ = %sDur$_{sub}$ + %s'%(slope, intercept), 
         color='green',
         #label = r'$\frac{d\nu}{dt}_{happy}$'
         )

#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Duration (ms)', fontsize=10)
plt.ylabel('Sub-Burst Drift Rate ' + r'$\frac{d\nu}{dt}$' + ' (MHz/ms)', fontsize=10)
#plt.legend()
plt.savefig('duration_sub_drift.pdf', dpi=600, bbox_inches='tight')
print('duration_sub_drift.png')