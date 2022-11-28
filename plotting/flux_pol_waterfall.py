import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import os
from scipy.signal import find_peaks

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), L (linear), SI (Stokes I), SV (Stokes V)

a = psrchive.Archive_load(sys.argv[1])
# polarisation type I,L,SI,V
p = sys.argv[2]

## DETERMINE PEAK FLUX AND INDEX FOR PLOT CENTRING
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


if sys.argv[2] == "I":
    ### ZOOMED POLARISATION
    c1 = a.clone()
    c1.remove_baseline()
    c1.tscrunch()
    c1.fscrunch(4)
    c1.pscrunch()
    data1 = c1.get_data()
    nsub, npol, nchan, nbin = data1.shape

    # on-pulse phase bin start and finish
    ps = int(p1*nbin)
    pf = int(p2*nbin)
    # zoom by phase on each side
    dp = 0.0485
    # on-pulse phase bin start and finish
    p3 = p1+dp
    p4 = p2-dp
    # on-pulse phase start and finish
    psz = int(p3*nbin)
    pfz = int(p4*nbin)

    # zoomed intensity
    Iz = data1[0,0,:,psz:pfz]


    ### POLARISATION
    c1.bscrunch(8)
    data2 = c1.get_data()
    nsub, npol, nchan, nbin = data2.shape

    # on-pulse phase bin start and finish
    ps = int(p1*nbin)
    pf = int(p2*nbin)

    # intensity
    I = data2[0,0,:,ps:pf]

else:
    ### ZOOMED POLARISATION
    c1 = a.clone()
    c1.remove_baseline()
    c1.tscrunch()
    c1.fscrunch(4)
    data1 = c1.get_data()
    nsub, npol, nchan, nbin = data1.shape

    # on-pulse phase bin start and finish
    ps = int(p1*nbin)
    pf = int(p2*nbin)
    # zoom by phase on each side
    dp = 0.0485
    # on-pulse phase bin start and finish
    p3 = p1+dp
    p4 = p2-dp
    # on-pulse phase start and finish
    psz = int(p3*nbin)
    pfz = int(p4*nbin)

    # zoomed polarisations
    SIz = data1[0,0,:,psz:pfz]
    Lz = np.sqrt(data1[0,1,:,psz:pfz]**2+data1[0,2,:,psz:pfz]**2)
    SVz = np.fabs(data1[0,3,:,psz:pfz])

    ### POLARISATION
    c1.bscrunch(8)
    data2 = c1.get_data()
    nsub, npol, nchan, nbin = data2.shape

    # on-pulse phase bin start and finish
    ps = int(p1*nbin)
    pf = int(p2*nbin)

    # polarisations
    SI = data2[0,0,:,ps:pf]
    L = np.sqrt(data2[0,1,:,ps:pf]**2+data2[0,2,:,ps:pf]**2)
    SV = np.fabs(data2[0,3,:,ps:pf])



### PLOT ZOOMED POLARISATION
ticks = np.round(np.linspace(p3,p4,num=11),4)
ticks_x = np.linspace(0,pfz-psz,num=11)

plt.figure(figsize=(15,15),dpi=300)
ax3 = plt.subplot(313)
plt.imshow(eval(p+'z'), cmap='Spectral', vmin=np.min(eval(p)), vmax=np.max(eval(p)), aspect='auto', origin='lower')
plt.xticks(ticks_x, ticks)
plt.xlabel('Phase')
plt.ylabel('Frequency channel')


### PLOT POLARISATION
ticks = np.round(np.linspace(p1,p2,num=11),4)
ticks_x = np.linspace(0,pf-ps,num=11)

ax1 = plt.subplot(312)
plt.imshow(eval(p), cmap='Spectral', vmin=np.min(eval(p)), vmax=0.3*np.max(eval(p)), aspect='auto', origin='lower')
plt.xticks(ticks_x, ticks)
plt.ylabel('Frequency channel')


### PLOT FLUX DENSITY
nsub, npol, nchan, nbin = data.shape

# on-pulse phase start and finish
ps = int(p1*nbin)
pf = int(p2*nbin)

ticks = np.round(np.linspace(p1,p2,num=11),4)
ticks_x = np.linspace(0,pf-ps,num=11)

ax2 = plt.subplot(311)
plt.plot(data[0,0,0,ps:pf], c='black')
ax2.set_xlim(0,pf-ps)
plt.xticks(ticks_x, ticks)
plt.ylabel('Flux Density')
plt.title('%s Polarisation %s'%(p,sys.argv[1]))

### SAVE FIGURE
plt.savefig('f_pol_waterfall_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
print('f_pol_waterfall_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
