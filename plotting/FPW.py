import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
import psrchive
import os
from scipy.signal import find_peaks

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

a = psrchive.Archive_load(sys.argv[1])
# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[2]
# neutron star period (ms)
period = 1/float(sys.argv[3])


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
p1 = np.round(peak_idx/nbin - 0.1, 4)
p2 = np.round(peak_idx/nbin + 0.1, 4)


if sys.argv[2] == "I":
    c1 = a.clone()
    c1.remove_baseline()
    c1.tscrunch()
    #c1.fscrunch(4)
    c1.pscrunch()
    data1 = c1.get_data()
    nsub, npol, nchan, nbin = data1.shape

    # on-pulse phase bin start and finish
    ps = int(np.round(p1*nbin))
    pf = int(np.round(p2*nbin))

    ### POLARISATION
    c1.bscrunch(8)
    data2 = c1.get_data()
    nsub, npol, nchan, nbin = data2.shape

    # on-pulse phase bin start and finish
    ps = int(np.round(p1*nbin))
    pf = int(np.round(p2*nbin))

    # intensity
    I = data2[0,0,:,ps:pf]

else:
    c1 = a.clone()
    c1.remove_baseline()
    c1.tscrunch()
    #c1.fscrunch(4)
    data1 = c1.get_data()
    nsub, npol, nchan, nbin = data1.shape

    ### POLARISATION
    c1.bscrunch(8)
    data2 = c1.get_data()
    nsub, npol, nchan, nbin = data2.shape

    # on-pulse phase bin start and finish
    ps = int(np.round(p1*nbin))
    pf = int(np.round(p2*nbin))

    # polarisations
    SI = data2[0,0,:,ps:pf]
    SQ = data2[0,1,:,ps:pf]
    SU = data2[0,2,:,ps:pf]
    L = np.sqrt(data2[0,1,:,ps:pf]**2+data2[0,2,:,ps:pf]**2)
    SV = data2[0,3,:,ps:pf]


#### PLOTTING ####
fig = plt.figure(figsize=(10,15),dpi=300)
g = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[1,7], hspace=0.)

# seconds per bin
bs = 1*period/nbin
nbin_zoom = np.shape(eval(p))[1]

### PLOT POLARISATION
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.linspace(704,4032, num=14).astype(int)
yticks_y = np.linspace(0,nchan-1, len(yticks))

#masked_data = np.ma.masked_values(eval(p), 0.)
#cmap = matplotlib.cm.get_cmap("Spectral").copy()
#cmap.set_bad(color='white')
#mask zapped channels colour
#if len(sys.argv)==5:
#    pol = np.ma.masked_values(eval(p), 0.)
#    cmap = matplotlib.cm.get_cmap("Spectral").copy()
#    cmap.set_bad(color='white')
#else:
#    pol = eval(p)

ax1 = fig.add_subplot(g[1])
ax1.imshow(eval(p), cmap="Spectral", vmin=np.min(eval(p)), vmax=0.3*np.max(eval(p)), aspect='auto', origin='lower')
plt.xticks(xticks_x, xticks, fontsize=10)
plt.yticks(yticks_y, yticks, fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Frequency (MHz)', fontsize=10)


### PLOT FLUX DENSITY
nsub, npol, nchan, nbin = data.shape

# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))

xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.round(np.linspace(0,peak_flux/1000, num=4)).astype(int)
yticks_y = np.linspace(0,peak_flux/1000,num=len(yticks))

ax2 = plt.subplot(g[0])
ax2.plot(data[0,0,0,ps:pf]/1000, c='black')
ax2.set_xlim(0,pf-ps)
ax2.set_ylim(-5,peak_flux/1000+5)
plt.xticks(xticks_x, xticks, fontsize=10)
ax2.axes.xaxis.set_ticklabels([])
plt.yticks(yticks_y, yticks, fontsize=10)
plt.ylabel('Flux Density (Jy)', fontsize=10)
plt.title('%s Polarisation %s'%(p,sys.argv[1]))

### SAVE FIGURE
plt.savefig('FPW_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
print('FPW_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
