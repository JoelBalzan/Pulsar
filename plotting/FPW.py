import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import psrchive
import os
from scipy import interpolate
from matplotlib.colorbar import Colorbar

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[2]

# Phase zoom factor
z = 0.0004

if sys.argv[2] == "I":
	a = psrchive.Archive_load(sys.argv[1])
	a.remove_baseline()
	a.tscrunch()
	a.pscrunch()
	data2 = a.get_data()
	nsub, npol, nchan, nbin = data2.shape

	# peak and index
	peak_idx = np.argmax(np.mean(data2[0,0,:,:], axis=0))

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

	# intensity
	I = data2[0,0,fs:ff,ps:pf]
    
else:
	a = psrchive.Archive_load(sys.argv[1])
	a.remove_baseline()
	a.tscrunch()
	#c1.bscrunch(8)
	data2 = a.get_data()
	nsub, npol, nchan, nbin = data2.shape

	# peak and index
	peak_idx = np.argmax(data2.mean(axis=(1,2))[0])

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

	# polarisations
	if p == "SI":
		SI = data2[0,0,fs:ff,ps:pf]
	if p == "SQ":
		SQ = data2[0,1,fs:ff,ps:pf]
	if p == "SU":
		SU = data2[0,2,fs:ff,ps:pf]
	if p == "L":
		L = np.sqrt(data2[0,1,fs:ff,ps:pf]**2+data2[0,2,fs:ff,ps:pf]**2)
	if p == "SV":
		SV = data2[0,3,fs:ff,ps:pf]

#### PLOTTING ####
fig = plt.figure(figsize=(10,15),dpi=300)
g = gridspec.GridSpec(ncols=2, nrows=2, height_ratios=[1,7], width_ratios=[14,1], hspace=0.,wspace=0.1)

# seconds per bin
period = a.integration_length()
bs = 1000*period/nbin
nbin_zoom = np.shape(eval(p))[1]

### PLOT DYNAMIC SPECTRUM
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.linspace(f1,f2, num=14).astype(int)
yticks_y = np.linspace(0,ff-fs-1, len(yticks))


ax1 = plt.subplot(g[1,0])
im = ax1.imshow(eval(p), cmap="viridis", vmin=np.min(eval(p)), 
	   vmax=0.7*np.max(eval(p)), aspect='auto', origin='lower', interpolation='kaiser')
ax1.set_xlim(0.0, pf-ps-1)
ax1.set_xticks(xticks_x)
ax1.set_xticklabels(xticks, fontsize=12)
ax1.set_xlabel('Time (ms)', fontsize=12)

ax1.set_ylabel('Frequency (MHz)', fontsize=12)
ax1.set_ylim(0, ff-fs-1)
ax1.set_yticks(yticks_y)
ax1.set_yticklabels(yticks, fontsize=12)

cbax = plt.subplot(g[1,1])
cb = Colorbar(ax=cbax, mappable=im, orientation='vertical')

### PLOT FLUX DENSITY
c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape

peak_flux = np.max(data[0,0,0,:])

# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))

xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.round(np.linspace(0,peak_flux/1000 - 1, num=4)).astype(int)
yticks_y = np.linspace(0,(peak_flux/1000)-1,num=len(yticks))

ax2 = plt.subplot(g[0,0])
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
plt.savefig('FPW_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]), bbox_inches='tight')
print('FPW_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
#plt.show()


