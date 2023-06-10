import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import psrchive
import os
from scipy.signal import find_peaks, peak_widths

# python pulse_drift.py <file> <Polarisation> <freq_window> <poly_degree>
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)


def cal_fwtm (freq, spec):
	peak = np.amax(spec)
	peak_idx = np.argmax(spec)
	tenth = 0.1*peak
	signs = np.sign(np.add(spec, -tenth))
	#print (peak, tenth, signs)
	zero_crossings = (signs[0:-1] != signs[1:])
	#zero_crossings = (signs[0:-2] != signs[1:-1])
	zero_crossings_i = np.where(zero_crossings)[0]
	
	#print (peak_idx)
	#print (zero_crossings, np.where(zero_crossings), zero_crossings_i, peak_idx)
	if len(zero_crossings_i) > 0:
		signs = np.sign(np.add(zero_crossings_i, -peak_idx))
		num = len(signs)
		#print (signs)
		if np.all(np.equal(signs, np.ones(num))):
			idx1 = 0
			idx2 = zero_crossings_i[0]
			#print ('here1')
		elif np.all(np.equal(signs, -np.ones(num))):
			idx1 = zero_crossings_i[-1]
			idx2 = -1
			#print ('here2')
		else:
			#print ('here3')
			#print (signs)
			#print (signs[0:-1], signs[1:])
			crossings2 = (signs[0:-1] != signs[1:])
			#crossings2 = (signs[0:-2] != signs[1:-1])
			#print (crossings2)
			crossings2_i = np.where(crossings2)[0]
			#print (crossings2_i)
		
			idx1 = zero_crossings_i[crossings2_i[0]]
			idx2 = zero_crossings_i[crossings2_i[0]+1]
	else:
		idx1 = 0
		idx2 = -1

	#print (idx1, idx2)
	freq1 = freq[idx1]
	freq2 = freq[idx2]
	#print (freq1, freq2)
	bw = np.fabs(freq1-freq2)

	#return freq1, freq2, np.fabs(freq1-freq2)
	return int(idx1), int(idx2), bw	



### DYNAMIC SPECTRA | PEAKS AND MINIMAS | WEIGHTED SPECTRA PEAKS ###
a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.pscrunch()
data2 = a.get_data()
nsub, npol, nchan, nbin = data2.shape
# peak and index
F = np.mean(data2[0,0,:,:], axis=0)/1000
peak_idx = np.argmax(F)
# peak width
width = np.round(8*peak_widths(F, np.array([peak_idx]), rel_height=0.5)[0]).astype(int)
if width%2 == 1:
	width = width - 1
# on-pulse phase start and finish
ps = int(np.round(peak_idx - width/2))
pf = int(np.round(peak_idx + width/2))
# dynamic spectra
P = data2[0,0,:,ps:pf]/1000

# PEAKS AND MINIMAS
h=0.3*np.amax(F[ps:pf])
peaks, _ = find_peaks(F[ps:pf], height=h)
# peak minimas
mins, _ = find_peaks(-F[ps:pf])
# associate peaks with minimas
if peaks[-1] > mins[-1]:
	peaks = peaks[:-1]
if peaks[0] < mins[0]:
	peaks = peaks[1:]
peak_mins = []
for i in peaks:
	for j in range(len(mins)):
		if mins[j] < i < mins[j+1]:
			mins_i = np.array([mins[j], mins[j+1]])
			peak_mins.append(mins_i)
peak_mins = np.array(peak_mins)
### SPECTRA OF PEAKS
spectra = data2[0,0,:,ps:pf]/1000
### DEFINE SPECTRA VARIABLES
S = []
for i in range(len(peaks)):
	S.append(np.mean(spectra[:,peak_mins[i][0]:peak_mins[i][1]], axis=1))
S = np.array(S)
### WEIGHTED SPECTRA CENTRE FREQUENCY (INDEX)
f_ch = np.array([a.get_first_Integration().get_centre_frequency(i) for i in range(nchan)])
bw = a.get_bandwidth()
cf = a.get_centre_frequency()
#lowest observed frequency
min_freq = cf-bw/2
# fscrunching factor
f_scr = bw/a.get_nchan()
spectra_centres = []
for i in range(len(peaks)):
	F_c = np.round((np.sum(np.multiply(S[i], range(nchan)))/np.sum(S[i])), 0).astype(int)
	spectra_centres.append(F_c)
spectra_centres = np.array(spectra_centres)



#### PLOTTING ####
A4x = 8.27
A4y = 11.69
fig = plt.figure(figsize=(1.5*A4x,1.5*A4y),dpi=300)
g = gridspec.GridSpec(ncols=2, nrows=2, height_ratios=[1,7], width_ratios=[7,1], hspace=0., wspace=0.)

# seconds per bin
period = a.integration_length()
bs = 1000*period/nbin
nbin_zoom = np.shape(P)[1]

### PLOT DYNAMIC SPECTRUM
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.linspace(704,4032, num=14).astype(int)
yticks_y = np.linspace(0,4032-704-1, len(yticks))


# use the variance to set the colour scale
var = np.var(P)/np.max(P)
vmin = 0.5*np.min(P)
vmax = 0.7*np.max(P)

ax1 = plt.subplot(g[1,0])
im = ax1.imshow(P, cmap="viridis", vmin=vmin, 
	   vmax=vmax, aspect='auto', origin='lower', interpolation='none')
ax1.scatter(peaks, spectra_centres, marker='.', c='k')
# plot best fit line
# Don't plot drift rate if there is only 1 peak detected
if (len(peaks) > 1):
	drift_rate_line = np.poly1d(np.polyfit(peaks, spectra_centres, 1))(np.unique(peaks))
	slope, intercept = np.polyfit(peaks, spectra_centres, 1)
	ax1.plot(np.unique(peaks), drift_rate_line, linestyle='-', c='k', 
		label=r'$\frac{d\nu}{dt}$ = %s MHz ms$^{-1}$'%(np.round(slope*((4032-704)/nchan)/(1000*period/nbin),2)))
	ax1.legend(loc='upper left')

ax1.set_xlim(0.0, pf-ps-1)
ax1.set_xticks(xticks_x)
ax1.set_xticklabels(xticks, fontsize=12)
ax1.set_xlabel('Time (ms)', fontsize=12)

ax1.set_ylabel('Frequency (MHz)', fontsize=12)
ax1.set_ylim(0, 4032-704-1)
ax1.set_yticks(yticks_y)
ax1.set_yticklabels(yticks, fontsize=12)
ax1.tick_params(axis='x', pad=10)


### PLOT PULSE PROFILE
peak_flux = np.max(F[ps:pf])

# on-pulse phase start and finish
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.round(np.linspace(0,peak_flux - 1, num=4)).astype(int)
yticks_y = np.linspace(0,(peak_flux)-1,num=len(yticks))

ax2 = plt.subplot(g[0,0])
#ax2.plot(np.arange(nbin), 0*np.arange(nbin), ls='--', color='gray')
ax2.plot(F[ps:pf], c='k', lw=1)
ax2.set_xlim(0,pf-ps-1)
ax2.set_ylim(-5,1.1*(peak_flux))
ax2.set_yticks([])
ax2.set_yticklabels([])
ax2.set_xticks(xticks_x)
ax2.set_xticklabels([])
ax2.set_ylabel('Flux Density (Jy)', fontsize=12, labelpad=20)
ax2.tick_params(axis="x", which='both', direction="in", pad=-10)
#plt.title('%s Polarisation %s'%(p,sys.argv[1].split('.')[0]))


### PLOT SPECTRA
spectrum = np.mean(P, axis=1)
ax3 = plt.subplot(g[1,1])

x = 0*np.arange(nchan)
ax3.plot(x, np.arange(nchan), ls='--', color='gray')
ax3.plot(spectrum, np.arange(nchan), ls='-', color='k', lw=1)

ax3.set_ylim(0.0, nchan-1)
#ax3.set_xlim(-3000,np.max(spectrum)+1000)
ax3.set_yticks(np.linspace(0,4032-704-1, num=14))
ax3.set_yticklabels([])
ax3.set_xticklabels([])
ax3.tick_params(axis="y", which='both', direction="in", pad=-22)


### SAVE FIGURE
plt.savefig('sub-pulse_drift_%s.pdf'%(sys.argv[1].split(os.extsep, 1)[0]), bbox_inches='tight')
print('sub-pulse_drift_%s.pdf'%(sys.argv[1].split(os.extsep, 1)[0]))
#plt.show()


