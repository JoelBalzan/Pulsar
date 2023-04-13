import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import psrchive
from matplotlib import gridspec
import os 
import sys
from scipy.signal import find_peaks, peak_widths


p = sys.argv[2]
if p == "I":
	a = psrchive.Archive_load(sys.argv[1])
	a.remove_baseline()
	a.tscrunch()
	a.pscrunch()
	data2 = a.get_data()
	nsub, npol, nchan, nbin = data2.shape

	# peak and index
	F = np.mean(data2[0,0,:,:], axis=0)
	peak_idx = np.argmax(F)

	# on-pulse phase start and finish
	width = np.round(4*peak_widths(F, np.array([peak_idx]), rel_height=0.5)[0]).astype(int)
	#p1 = np.round(peak_idx/nbin - width/2, 4)
	#p2 = np.round(peak_idx/nbin + width/2, 4)

	# on-pulse phase bin start and finish
	ps = int(np.round(peak_idx - width/2))
	pf = int(np.round(peak_idx + width/2))

	### FREQ ZOOM
	f_scr = (4032-704)/a.get_nchan()

	f1 = float(sys.argv[3])
	f2 = float(sys.argv[4])
	fs = int((f1-704)/f_scr)
	ff = int((f2-704)/f_scr)

	# intensity
	P = data2[0,0,fs:ff,ps:pf]

else:
	a = psrchive.Archive_load(sys.argv[1])
	a.remove_baseline()
	a.tscrunch()
	#a.bscrunch(2)
	data2 = a.get_data()
	nsub, npol, nchan, nbin = data2.shape

	# peak and index
	F = data2.mean(axis=(1,2))[0]
	peak_idx = np.argmax(F)

	# on-pulse phase start and finish
	width = np.round(4*peak_widths(F[peak_idx], peak_idx, rel_height=0.5)[0]).astype(int)
	#p1 = np.round(peak_idx/nbin - width/2, 4)
	#p2 = np.round(peak_idx/nbin + width/2, 4)

	# on-pulse phase bin start and finish
	ps = int(np.round(peak_idx - width/2))
	pf = int(np.round(peak_idx + width/2))

	### FREQ ZOOM
	f_scr = (4032-704)/a.get_nchan()

	f1 = float(sys.argv[3])
	f2 = float(sys.argv[4])
	fs = int((f1-704)/f_scr)
	ff = int((f2-704)/f_scr)

	# polarisations
	if p == "SI":
		P = data2[0,0,fs:ff,ps:pf]
	if p == "SQ":
		P = data2[0,1,fs:ff,ps:pf]
	if p == "SU":
		P = data2[0,2,fs:ff,ps:pf]
	if p == "L":
		P = np.sqrt(data2[0,1,fs:ff,ps:pf]**2+data2[0,2,fs:ff,ps:pf]**2)
	if p == "SV":
		P = data2[0,3,fs:ff,ps:pf]

# 1D autocorrelation of each frequency channel
nchan, nbin = P.shape

freq_corr = []
for i in range(nchan):
	freq_corr.append(signal.correlate(P[i,:], P[i,:], mode='full', method='direct'))

# 1D autocorrelation of each phase bin
phase_corr = []
for i in range(nbin):
	phase_corr.append(signal.correlate(P[:,i], P[:,i], mode='full', method='direct'))
phase_corr = np.rot90(phase_corr)

# 2D autocorrelation
corr_2D = signal.correlate2d(P, P, mode='full', boundary='fill', fillvalue=0)

### PLOTTING ###
A4x, A4y = 8.27, 11.69
fig = plt.figure(figsize=(A4x, A4x), dpi=300)
g = gridspec.GridSpec(ncols=3, nrows=3, hspace=0., wspace=0.)

# plot dynamic spectrum
ax_2_0 = fig.add_subplot(g[2,0])
ax_2_0.imshow(P, cmap='Greys', aspect='auto', origin='lower', interpolation='none')

# plot frequency autocorrelation
ax_2_1 = fig.add_subplot(g[2,1])
ax_2_1.imshow(freq_corr, cmap='Reds', aspect='auto', origin='lower', interpolation='none')

# plot phase autocorrelation
ax_1_0 = fig.add_subplot(g[1,0])
ax_1_0.imshow(phase_corr, cmap='Blues', aspect='auto', origin='lower', interpolation='none')

# plot 2D autocorrelation
ax_1_1 = fig.add_subplot(g[1,1])
ax_1_1.imshow(corr_2D, cmap='Purples', aspect='auto', origin='lower', interpolation='none')

plt.savefig(sys.argv[1]+'.png', dpi=300, bbox_inches='tight')
