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
#p1 = np.round(peak_idx/nbin - 0.1, 4)
#p2 = np.round(peak_idx/nbin + 0.1, 4)
p1 = np.round(peak_idx/nbin - 0.0008, 4)
p2 = np.round(peak_idx/nbin + 0.0008, 4)


if sys.argv[2] == "I":
	c1 = a.clone()
	c1.remove_baseline()
	c1.tscrunch()
	c1.pscrunch()
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
bs = 1000*period/nbin
nbin_zoom = np.shape(eval(p))[1]

### PLOT DYNAMIC SPECTRUM
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.linspace(704,4032, num=14).astype(int)
yticks_y = np.linspace(0,nchan, len(yticks))

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
ax1.set_xlim(0.0, pf-ps-1)
ax1.set_xticks(xticks_x)
ax1.set_xticklabels(xticks, fontsize=12)
ax1.set_xlabel('Time (ms)', fontsize=12)

ax1.set_ylabel('Frequency (MHz)', fontsize=12)
ax1.set_ylim(0, nchan-1)
ax1.set_yticks(yticks_y)
ax1.set_yticklabels(yticks, fontsize=12)



### PLOT FLUX DENSITY
nsub, npol, nchan, nbin = data.shape

# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))

xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.round(np.linspace(0,peak_flux/1000 - 1, num=4)).astype(int)
yticks_y = np.linspace(0,(peak_flux/1000)-1,num=len(yticks))

ax2 = plt.subplot(g[0])
ax2.plot(data[0,0,0,ps:pf]/1000, c='black')
ax2.plot(np.arange(nbin), 0*np.arange(nbin), ls='--', color='k')
ax2.set_xlim(0,pf-ps-1)
ax2.set_ylim(-5,peak_flux/1000+5)
ax2.set_xticks(xticks_x)
ax2.set_xticklabels([])
ax2.set_ylabel('Flux Density (Jy)', fontsize=12)
ax2.tick_params(axis="x", which='both', direction="in", pad=-22)
plt.title('%s Polarisation %s'%(p,sys.argv[1].split('.')[0]))

### SAVE FIGURE
plt.savefig('FPW_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
print('FPW_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
