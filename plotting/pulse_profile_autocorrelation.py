import os
import sys

import lmfit
import matplotlib.pyplot as plt
import numpy as np
import psrchive
from scipy import signal
from scipy.signal import peak_widths

if os.path.isfile('corr_1D_%s.npy'%sys.argv[1].split(os.extsep, 1)[0]):
    corr_1D = np.load('corr_1D_%s.npy'%sys.argv[1].split(os.extsep, 1)[0])
else:
    a = psrchive.Archive_load(sys.argv[1])
    a.remove_baseline()
    a.tscrunch()
    a.fscrunch()
    a.pscrunch()
    data = a.get_data()
    nsub, npol, nchan, nbin = data.shape

    corr_1D = signal.correlate(data[0,0,0,:], data[0,0,0,:], mode='full', method='direct')
    np.save('corr_1D_%s.npy'%sys.argv[1].split(os.extsep, 1)[0], corr_1D)

A4x, A4y = 8.27, 11.69
fig = plt.figure(figsize=(A4y, A4x), dpi=600)

plt.plot(corr_1D[int(len(corr_1D)/2):int(len(corr_1D)/2)+5000], c='k', lw=0.5)

plt.savefig(sys.argv[0].split(os.extsep, 1)[0]+'_%s_'%sys.argv[1]+sys.argv[1].split(os.extsep, 1)[0]+'.png', dpi=600, bbox_inches='tight')
print(sys.argv[0].split(os.extsep, 1)[0]+'_%s_'%sys.argv[1]+sys.argv[1].split(os.extsep, 1)[0]+'.png')