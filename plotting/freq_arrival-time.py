import psrchive
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks, peak_widths
import sys

a = psrchive.Archive_load(sys.argv[1])

### DETERMINE PEAK FLUX AND INDEX FOR PLOT CENTRING
c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape

# peak and index
flux = data[0,0,0,:]
peaks, _ = find_peaks(flux)
peak_flux = np.sort(flux[peaks])[-1]
peak_idx = np.where(flux==peak_flux)[0][0]

# on-pulse phase start and finish
#p1 = np.round(peak_idx/nbin - 0.0008, 4)
#p2 = np.round(peak_idx/nbin + 0.0008, 4)
p1 = np.round(peak_idx/nbin - 0.001, 4)
p2 = np.round(peak_idx/nbin + 0.001, 4)

### PEAK INDICES  
# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))


a.remove_baseline()
a.tscrunch()
a.fscrunch(16)
a.pscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

# peak flux in each freq channel
peak_flux = np.array([np.argmax(data[0,0,i,ps:pf]) for i in range(nchan)]).astype(float)
peak_flux[peak_flux == 0.] = np.nan
f_ch = np.array([a.get_first_Integration().get_centre_frequency(i) for i in range(nchan)])

### CURVE FITTING ###


### PLOTTING ###
xticks = np.linspace(704,4032, num=14).astype(int)
xticks_x = np.linspace(0,nchan-1, len(xticks))

fig = plt.figure(figsize=(15, 10), dpi=300) 
plt.plot(peak_flux, color='k', lw=1, label='Arrival time')
#plt.plot(1/f_ch**2, color='b', lw=1, ls='--', label=r'$\frac{1}{\nu^{2}}$')
#plt.plot(1/f_ch**4, color='r', lw=1, ls='--', label=r'$\frac{1}{\nu^{4}}$')
plt.xlabel('Frequency (MHz)')
plt.xlim(0,nchan-1)
plt.xticks(xticks_x, xticks)
plt.ylabel('')
plt.legend()
plt.title('%s'%sys.argv[1].split('.')[0], fontsize=12)

plt.savefig("freq_arrival-time_%s.pdf"%(sys.argv[1].split('.')[0]))
print("freq_arrival-time_%s.pdf"%(sys.argv[1].split('.')[0]))