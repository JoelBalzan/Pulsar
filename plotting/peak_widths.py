import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
from scipy.signal import find_peaks, peak_widths
import os

# python peak_widths.py <file> 


a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.fscrunch()
a.pscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

#spin period of object
period = a.integration_length()


# peak and index
flux = data[0,0,0,:]
peaks, _ = find_peaks(flux)
peak_flux = np.sort(flux[peaks])[-1]
peak_idx = np.where(flux==peak_flux)[0][0]


# on-pulse phase start and finish
p1 = np.round(peak_idx/nbin - 0.05, 4)
p2 = np.round(peak_idx/nbin + 0.05, 4)
ps = int(p1*nbin)
pf = int(p2*nbin)

# flux in Jy
flux = data[0,0,0,ps:pf]/1000
peaks, _ = find_peaks(flux)

# highest 10 fluxes
fluxes = np.sort(flux[peaks])[::-1][0:10]

# index of highest 10 fluxes
fluxes_i = []
for i in fluxes:
    idx = np.where(flux==i)[0][0]
    fluxes_i.append(idx)
fluxes_i = np.array(fluxes_i)

# width
w = peak_widths(flux, fluxes_i, rel_height=0.98)


## FLUX DENSITY PLOT WITH PEAKS AND WIDTHS MARKED
ticks = np.round(np.arange(p1,p2+0.01, step=0.01), 2)
ticks_x = np.linspace(0,pf-ps,num=len(ticks))

plt.figure(figsize=(15,10),dpi=300)
ax = plt.subplot(111)
ax.set_xlim(0,pf-ps)
plt.plot(flux, label='Flux Density')
plt.plot(fluxes_i, flux[fluxes_i], 'x', label='Peaks')
plt.hlines(*w[1:], color="C2")
plt.xticks(ticks_x, ticks)
plt.xlabel('Phase')
plt.ylabel('Flux Density')
plt.legend()

plt.savefig("peaks_widths_%s.pdf"%(sys.argv[1].split(os.extsep, 1)[0]))

