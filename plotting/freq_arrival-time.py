import psrchive
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks, peak_widths
import sys

a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.bscrunch(2)
a.pscrunch()
a.centre()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

# peak and index
peak_idx = np.argmax(np.mean(data[0,0,:,:], axis=0))

# Phase zoom factor
z = 0.0004
# on-pulse phase start and finish
p1 = np.round(peak_idx/nbin - z, 4)
p2 = np.round(peak_idx/nbin + z, 4)

# on-pulse phase bin start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))

### FREQ ZOOM
f_scr = (4032-704)/a.get_nchan()

f1 = float(sys.argv[3])
f2 = float(sys.argv[4])
fs = int((f1-704)/f_scr)
ff = int((f2-704)/f_scr)


# peak flux in each freq channel
peak_flux = np.array([np.argmax(data[0,0,i,ps:pf]) for i in np.arange(fs,ff)]).astype(float)
peak_flux[peak_flux == 0.] = np.nan
#f_ch = np.array([a.get_first_Integration().get_centre_frequency(i) for i in range(nchan)])

### CURVE FITTING ###


### PLOTTING ###
xticks = np.linspace(f1,f2, num=14).astype(int)
xticks_x = np.linspace(0,ff-fs-1, len(xticks))

fig = plt.figure(figsize=(15, 10), dpi=300) 
x = np.arange(0,ff-fs)
plt.scatter(x, peak_flux, color='k', s = 1)
#plt.plot(1/f_ch**2, color='b', lw=1, ls='--', label=r'$\frac{1}{\nu^{2}}$')
#plt.plot(1/f_ch**4, color='r', lw=1, ls='--', label=r'$\frac{1}{\nu^{4}}$')
plt.xlabel('Frequency (MHz)')
plt.xticks(xticks_x, xticks)
plt.yticks([])
plt.ylabel('Arrival time')
plt.title('%s'%sys.argv[1].split('.')[0], fontsize=12)

plt.savefig("freq_arrival-time_%s.pdf"%(sys.argv[1].split('.')[0]))
print("freq_arrival-time_%s.pdf"%(sys.argv[1].split('.')[0]))