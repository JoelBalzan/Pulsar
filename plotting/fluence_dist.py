import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob
from scipy.signal import find_peaks, peak_widths

# python energy_dist.py period


# top 10 fluxes
fluxes = np.empty(0)
# indexes of top 10 fluxes
fluxes_i = np.empty(0)
#widths in ms
widths = np.empty(0)
#spin period of object
period = 1/float(sys.argv[1])

counter = 0

for ar in glob.glob("*.rescaled"):
    #a = psrchive.Archive_load(sys.argv[1])
    a = psrchive.Archive_load(ar)
    a.remove_baseline()
    a.tscrunch()
    a.fscrunch()
    a.pscrunch()
    data = a.get_data()
    nsub, npol, nchan, nbin = data.shape

    
    # peak and index
    flux = data[0,0,0,:]
    peaks, _ = find_peaks(flux)
    peak_flux = max(flux[peaks])
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
    w = peak_widths(flux, fluxes_i[-10:], rel_height=0.95)[0]*bs
    widths = np.concatenate((widths,w), axis=None)

    counter = counter + 1
    print("%s/%s"%(counter,len(glob.glob("*.rescaled")))," files completed", end='\r')
    
fluence = fluxes*widths

## PLOT HISTOGRAM
hist, bins, _ = plt.hist(fluence, bins=30, histtype='step')
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

plt.figure(figsize=(15,15),dpi=300)
ax = plt.subplot(111)
plt.hist(fluence, bins=logbins, histtype='step')
ax.set_xscale("log")
ax.set_yscale("log") 
plt.xlabel('log$_{10}$(Fluence) (Jy ms)')
plt.ylabel('log$_{10}$(Count)')
plt.title('PX500 38329')

plt.savefig("e_dist_PX500_loglog.png")


