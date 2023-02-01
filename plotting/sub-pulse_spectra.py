import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import psrchive
import os
from scipy.signal import find_peaks, peak_widths
from scipy.signal import savgol_filter



a = psrchive.Archive_load(sys.argv[1])


### DETERMINE PEAK FLUX AND INDEX FOR PLOT CENTRING
c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape

# neutron star period (ms)
period = a.integration_length()


# peak and index
flux = data[0,0,0,:]
peaks, _ = find_peaks(flux)
peak_flux = np.sort(flux[peaks])[-1]
peak_idx = np.where(flux==peak_flux)[0][0]

# on-pulse phase start and finish
p1 = np.round(peak_idx/nbin - 0.001, 4)
p2 = np.round(peak_idx/nbin + 0.001, 4)


### PEAK INDICES  
# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))

flux = data[0,0,0,ps:pf]/1000
peaks, _ = find_peaks(flux)

# highest fluxes
fluxes = np.sort(flux[peaks])[::-1][0:6]

# index of highest fluxes
fluxes_i = []
for i in fluxes:
    idx = np.where(flux==i)[0][0]
    fluxes_i.append(idx)
fluxes_i = np.array(fluxes_i)

# peak widths
widths = peak_widths(flux, fluxes_i, rel_height=0.4)

# peak minimas
mins, _ = find_peaks(-flux)
# associate peaks with minimas
peak_mins = []
for i in fluxes_i:
    for j in range(len(mins)):
        if mins[j] < i < mins[j+1]:
            mins_i = np.array([[mins[j], mins[j+1]]])[0]
            peak_mins.append(mins_i)
peak_mins = np.array(peak_mins)


fluxes_i, fluxes = zip(*sorted(zip(fluxes_i, fluxes)))

### SPECTRA OF PEAKS
c1 = a.clone()
c1.remove_baseline()
c1.tscrunch()
c1.pscrunch()
data1 = c1.get_data()
nsub, npol, nchan, nbin = data1.shape
spectra = data1[0,0,:,ps:pf]/1000

dict={}
for i in range(len(fluxes)):
    key = str("S"+str(i))
    spec = np.mean(spectra[:,peak_mins[i][0]:peak_mins[i][1]], axis=1)
    dict[key] = spec.tolist()
for key,value in dict.items():
    exec(f'{key}={value}')


#### PLOTTING ####
fig = plt.figure(figsize=(15,15),dpi=300)
g = gridspec.GridSpec(ncols=1, nrows=len(fluxes), hspace=0)

freq_window = float(sys.argv[2])
freq_poly = int(sys.argv[3])
window = freq_window*nchan
window = int(np.ceil(window) // 2*2 +1)

for i in range(len(fluxes)):
    ax = fig.add_subplot(g[i])
    ax.plot(eval("S"+str(i)), c='k', label='peak %s'%i)
    fit_spec = savgol_filter(eval("S"+str(i)), window, freq_poly)
    ax.plot(np.arange(nchan), fit_spec, ls='--', color='r', label='Deg = %s, Window = %s'%(freq_poly, freq_window))
    ax.set_xlim(0,nchan-1)
    ax.legend(loc='upper left')
ax.set_xticks(np.linspace(0,nchan-1, 14))
ax.set_xticklabels(np.linspace(704,4032, num=14).astype(int))
ax.set_xlabel('Frequency (MHz)')


### SAVE FIGURE
plt.savefig("sub-pulse_spectra_%s.pdf"%(sys.argv[1].split('.')[0]), bbox_inches='tight')
print("sub-pulse_spectra_%s.pdf"%(sys.argv[1].split('.')[0]))


