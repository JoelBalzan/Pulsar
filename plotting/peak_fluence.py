import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob
from scipy.signal import find_peaks, peak_widths

# python peak_fluence.py <file> <period>


a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.fscrunch()
a.pscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

# peak and index
flux = data[0,0,0,:]
peaks, _ = find_peaks(flux)
peak_flux = np.sort(flux[peaks])[-1]
peak_idx = np.where(flux==peak_flux)[0][0]
# on-pulse phase start and finish
p1 = np.round(peak_idx/nbin - 0.1, 4)
p2 = np.round(peak_idx/nbin + 0.1, 4)
ps = int(p1*nbin)
pf = int(p2*nbin)
flux = data[0,0,0,ps:pf]/1000
peaks, _ = find_peaks(flux)

# highest fluxes
flux1 = np.sort(flux[peaks])[::-1][0]

# index of highest 10 fluxes
idx = np.array([np.where(flux==flux1)[0][0]])

# milliseconds per bin
period = a.integration_length()
bs = 1000*period/nbin
# width in milliseconds
width = peak_widths(flux, idx, rel_height=0.98)[0]*bs


fluence = flux1*width[0]
print("%s = %s Jy ms"%(sys.argv[1],fluence))
