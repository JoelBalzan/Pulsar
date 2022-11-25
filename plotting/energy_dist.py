import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob
import os
import pandas as pd
from scipy.signal import find_peaks, peak_widths

# top 10 fluxes
fluxes = np.empty(0)
# indexes of top 10 fluxes
fluxes_i = np.empty(0)
#widths in ms
widths = np.empty(0)
#spin period of object
period = 1/float(sys.argv[3])


for ar in glob.glob("*.rescaled"):
    #a = psrchive.Archive_load(sys.argv[1])
    a = psrchive.Archive_load(ar)
    a.remove_baseline()
    a.tscrunch()
    a.fscrunch()
    a.pscrunch()
    data = a.get_data()
    nsub, npol, nchan, nbin = data.shape

    p1 = float(sys.argv[1])
    p2 = float(sys.argv[2])
    ps = int(p1*nbin)
    pf = int(p2*nbin)


    flux = data[0,0,0,ps:pf]
    peaks, _ = find_peaks(flux)
    # highest 10 fluxes
    flux10 = np.sort(flux[peaks])[::-1][0:10]
    fluxes = np.concatenate((fluxes,flux10), axis=None)


    # index of highest 10 fluxes
    idxs = []
    for i in fluxes[-10:]:
        idx = np.where(flux==i)[0][0]
        idxs.append(idx)
    idxs = np.array(idxs)
    fluxes_i = np.concatenate((fluxes_i,idxs), axis=None).astype(int)


    # milliseconds per bin
    bs = 1000*period/nbin
    # width in milliseconds
    w = peak_widths(flux, fluxes_i[-10:], rel_height=0.9)[0]*bs
    widths = np.concatenate((widths,w), axis=None)


fluence = fluxes*widths


## PLOT HISTOGRAM
plt.figure(figsize=(15,15),dpi=300)
ax = plt.subplot(111)
plt.hist(fluence, bins=50, histtype='step')
plt.xlabel('Fluence (Jy ms)')
plt.ylabel('Count')
plt.title('P970')

plt.savefig("e_dist_P970.png")



'''
## FLUX DENSITY PLOT WITH PEAKS AND WIDTHS MARKED
ticks = np.round(np.arange(p1,p2+0.01, step=0.01), 2)
ticks_x = np.linspace(0,pf-ps,num=len(ticks))

plt.figure(figsize=(15,10),dpi=300)
plt.plot(flux, label='Flux Density')
plt.plot(flux10_i, flux[flux10_i], 'x', label='Peaks')
plt.hlines(*widths[1:], colors="C2", label='widths')
plt.xticks(ticks_x, ticks)
plt.xlabel('Phase')
plt.ylabel('Flux Density')
plt.legend()

plt.savefig("flux_peaks_widths_%s.png"%(sys.argv[1].split(os.extsep, 1)[0]))
'''
