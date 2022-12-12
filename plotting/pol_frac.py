import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
from scipy.signal import find_peaks

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
p1 = np.round(peak_idx/nbin - 0.05, 4)
p2 = np.round(peak_idx/nbin + 0.05, 4)


### POLARISATION FRACTION
c1 = a.clone()
c1.remove_baseline()
c1.tscrunch()
c1.fscrunch()
data1 = c1.get_data()
nsub, npol, nchan, nbin = data1.shape

ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))

I = data1[0,0,0,ps:pf]
I[I==0] = np.nan
Q = data1[0,1,0,ps:pf]
U = data1[0,2,0,ps:pf]
V = data1[0,3,0,ps:pf]

L = np.sqrt(np.float64(Q**2+U**2))/np.fabs(I)
C = np.fabs(np.float64(V)/I)


xticks = np.round(np.linspace(p1,p2,num=11),4)
xticks_x = np.linspace(0,pf-ps,num=len(xticks))

plt.figure(figsize=(15,10),dpi=300)
ax = plt.subplot(111)
#ax.set_ylim(0,1)
plt.plot(C, c='b', label='Circular')
plt.plot(L, c='r', label='Linear')
plt.xticks(xticks_x, xticks)
plt.xlabel('Phase')
plt.ylabel('Polarisation Fraction')
plt.legend()

plt.savefig("pol_frac.png")
