import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import psrchive
import os
from scipy import interpolate
from matplotlib.colorbar import Colorbar
from scipy.signal import find_peaks, peak_widths

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[2]


if sys.argv[2] == "I":
	a = psrchive.Archive_load(sys.argv[1])
	a.remove_baseline()
	a.tscrunch()
	a.pscrunch()
	data = a.get_data()
	nsub, npol, nchan, nbin = data.shape
	# peak and index
	pulse_profile = np.mean(data[0,0,:,:], axis=0)
	peak_idx = np.argmax(pulse_profile)
    # width of peaks for setting imshow widths
	width = peak_widths(pulse_profile, np.array([peak_idx]), rel_height=0.5)
	w = np.round(8*width[0]).astype(int)
	# on-pulse phase bin start and finish
	ps = int(np.round(peak_idx - w/2))
	pf = int(np.round(peak_idx + w/2))
	### FREQ ZOOM
	f1 = float(sys.argv[3])
	f2 = float(sys.argv[4])
	bw = a.get_bandwidth()
	cf = a.get_centre_frequency()
	#lowest observed frequency
	min_freq = cf-bw/2
	# fscrunching factor
	f_scr = bw/a.get_nchan()

	fs = int((f1-min_freq)/f_scr)
	ff = int((f2-min_freq)/f_scr)
	# intensity
	I = data[0,0,fs:ff,ps:pf]
	pulse_profile = pulse_profile[ps:pf]
else:
	a = psrchive.Archive_load(sys.argv[1])
	a.remove_baseline()
	a.tscrunch()
	a.centre()
	#a.bscrunch(2)
	data = a.get_data()
	nsub, npol, nchan, nbin = data.shape
	# peak and index
	pulse_profile = data.mean(axis=(1,2))[0]
	peak_idx = np.argmax(data.mean(axis=(1,2))[0])
    # width of peaks for setting imshow widths
	width = peak_widths(pulse_profile, np.array([peak_idx]), rel_height=0.8)
	w = np.round(3*width[0]).astype(int)
	# on-pulse phase bin start and finish
	ps = int(peak_idx - w)
	pf = int(peak_idx + w + 1)
	### FREQ ZOOM
	f_scr = (4032-704)/a.get_nchan()
	f1 = float(sys.argv[3])
	f2 = float(sys.argv[4])
	fs = int((f1-704)/f_scr)
	ff = int((f2-704)/f_scr)
	# polarisations
	if p == "SI":
		SI = data[0,0,fs:ff,ps:pf]
	if p == "SQ":
		SQ = data[0,1,fs:ff,ps:pf]
	if p == "SU":
		SU = data[0,2,fs:ff,ps:pf]
	if p == "L":
		L = np.sqrt(data[0,1,fs:ff,ps:pf]**2+data[0,2,fs:ff,ps:pf]**2)
	if p == "SV":
		SV = data[0,3,fs:ff,ps:pf]
	pulse_profile = pulse_profile[ps:pf]
#### PLOTTING ####
A4x = 8.27
A4y = 11.69
fig = plt.figure(figsize=(A4x,A4y),dpi=300)
g = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[1,7], hspace=0.)
# seconds per bin
period = a.integration_length()
bs = 1000*period/nbin
nbin_zoom = np.shape(eval(p))[1]
### PLOT DYNAMIC SPECTRUM
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.linspace(f1,f2, num=14).astype(int)
yticks_y = np.linspace(0,ff-fs-1, len(yticks))
# use the variance to set the colour scale
var = np.var(eval(p))/np.max(eval(p))
vmin = 0.5*np.min(eval(p))
vmax = 0.7*np.max(eval(p))
yticks = np.linspace(f1,f2, num=14).astype(int)
yticks_y = np.linspace(0,ff-fs-1, len(yticks))
ax1 = plt.subplot(g[1,0])
im = ax1.imshow(eval(p), cmap="viridis", 
	 vmin=vmin, 
	 vmax=vmax, 
	 aspect='auto', origin='lower', interpolation='none')
ax1.set_xlim(0.0, pf-ps-1)
ax1.set_xticks(xticks_x)
ax1.set_xticklabels(xticks, fontsize=10)
ax1.set_xlabel('Time (ms)', fontsize=10)
ax1.set_ylabel('Frequency (MHz)', fontsize=10)
ax1.set_ylim(0, ff-fs-1)
ax1.set_yticks(yticks_y)
ax1.set_yticklabels(yticks, fontsize=10)
ax1.tick_params(axis='x', pad=10)
#cbax = plt.subplot(g[1,1])
#cb = Colorbar(ax=cbax, mappable=im, orientation='vertical')
### PLOT PULSE PROFILE
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
ax2 = plt.subplot(g[0,0])
#ax2.plot(np.arange(nbin), 0*np.arange(nbin), ls='--', color='gray')
ax2.plot(pulse_profile, c='k', lw=1)
ax2.set_xlim(0,pf-ps-1)
ax2.set_yticks([])
ax2.set_yticklabels([])
ax2.set_xticks(xticks_x)
ax2.set_xticklabels([])
#ax2.set_ylabel('Flux Density (Jy)', fontsize=12, labelpad=20)
ax2.tick_params(axis="x", which='both', direction="in", pad=-10)
#plt.title('%s Polarisation %s'%(p,sys.argv[1].split('.')[0]))
### SAVE FIGURE
plt.savefig('FPW_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]), bbox_inches='tight')
print('FPW_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
#plt.show()


