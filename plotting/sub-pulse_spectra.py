import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import psrchive
import os
from scipy.signal import find_peaks, peak_widths


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


#### PLOTTING ####
fig = plt.figure(figsize=(10, 12), dpi=300) 

fluxes_i, freq_i = zip(*sorted(zip(fluxes_i, fluxes)))

### PLOT FLUX DENSITY
nsub, npol, nchan, nbin = data.shape

# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))



c1 = a.clone()
c1.remove_baseline()
c1.tscrunch()
c1.pscrunch()
data1 = c1.get_data()
nsub, npol, nchan, nbin = data1.shape
spectra = data1[0,0,:,ps:pf]/1000
print(spectra[:,28:34].shape)

dict={}
for i in range(len(fluxes)):
    key = str("S"+str(i))
    #dict[key] = spectra[:,peak_mins[i][0]:peak_mins[i][1]].tolist()
    dict[key] = 0
for key,value in dict.items():
    exec(f'{key}={value}')
print(len(S2))







#ax0.set_title('%s Polarisation %s'%(p,sys.argv[1].split('.')[0]), fontsize=12)




### SAVE FIGURE
#plt.savefig("pulse_drift_%s.pdf"%(sys.argv[1].split('.')[0]))
#print("pulse_drift_%s.pdf"%(sys.argv[1].split('.')[0]))


