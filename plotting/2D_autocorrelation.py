import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import psrchive
from matplotlib import gridspec
import os 
import sys
from scipy.signal import peak_widths


p = sys.argv[2]
if p == "I":
	a = psrchive.Archive_load(sys.argv[1])
	a.remove_baseline()
	a.tscrunch()
	a.pscrunch()
	data2 = a.get_data()
	nsub, npol, nchan, nbin = data2.shape

	# ms per bin
	period = a.integration_length()
	mspb = 1000*period/nbin

	# peak and index
	F = np.mean(data2[0,0,:,:], axis=0)
	peak_idx = np.argmax(F)

	# on-pulse phase start and finish
	width = np.round(6*peak_widths(F, np.array([peak_idx]), rel_height=0.5)[0]).astype(int)
	#p1 = np.round(peak_idx/nbin - width/2, 4)
	#p2 = np.round(peak_idx/nbin + width/2, 4)

	# on-pulse phase bin start and finish
	ps = int(np.round(peak_idx - width/2))
	pf = int(np.round(peak_idx + width/2))

	### FREQ ZOOM
	bw = a.get_bandwidth()
	cf = a.get_centre_frequency()
	#lowest observed frequency
	min_freq = cf-bw/2
	f_scr = bw/a.get_nchan()

	f1 = int(sys.argv[3])
	f2 = int(sys.argv[4])
	fs = int((f1-min_freq)/f_scr)
	ff = int((f2-min_freq)/f_scr)

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
	width = np.round(6*peak_widths(F[peak_idx], peak_idx, rel_height=0.5)[0]).astype(int)
	#p1 = np.round(peak_idx/nbin - width/2, 4)
	#p2 = np.round(peak_idx/nbin + width/2, 4)

	# on-pulse phase bin start and finish
	ps = int(np.round(peak_idx - width/2))
	pf = int(np.round(peak_idx + width/2))

	### FREQ ZOOM
	bw = a.get_bandwidth()
	cf = a.get_centre_frequency()
	f_scr = bw/a.get_nchan()

	f1 = int(sys.argv[3])
	f2 = int(sys.argv[4])
	fs = int((f1-(cf-bw/2))/f_scr)
	ff = int((f2-(cf-bw/2))/f_scr)

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
nchan, nbin = P.shape

## 1D auto-correlation of each frequency channel
freq_corr = []
for i in range(nchan):
	freq_corr.append(signal.correlate(P[i,:], P[i,:], mode='full', method='direct'))
# summed frequency auto-correlation
sum_freq_corr_freq = np.sum(freq_corr, axis=1)
sum_freq_corr_time = np.sum(freq_corr, axis=0)

## 1D auto-correlation of each phase bin
phase_corr = []
for i in range(nbin):
	phase_corr.append(signal.correlate(P[:,i], P[:,i], mode='full', method='direct'))
phase_corr = np.rot90(phase_corr)
# summed phase auto-correlation
sum_phase_corr_freq = np.sum(phase_corr, axis=1)
sum_phase_corr_time = np.sum(phase_corr, axis=0)

## 2D auto-correlation
corr_2D = signal.correlate2d(P, P, mode='full', boundary='fill', fillvalue=0)
sum_corr_2D_freq = np.sum(corr_2D, axis=1)
# summed phase auto-correlation
sum_corr_2D_freq[sum_corr_2D_freq==0] = np.nan
sum_corr_2D_time = np.sum(corr_2D, axis=0)


### PLOTTING ###
A4x, A4y = 8.27, 11.69
fig = plt.figure(figsize=(A4x, A4x), dpi=300)
g = gridspec.GridSpec(ncols=3, nrows=3, hspace=0., wspace=0., 
		      height_ratios=[0.5,1,1], width_ratios=[1,1,0.5])

# plot dynamic spectrum
ms_tick = nbin*mspb
dy_spec_xticks = np.round(np.linspace(0,ms_tick,num=5),2)
dy_spec_xticks_x = np.linspace(0,pf-ps-1,num=len(dy_spec_xticks))
dy_spec_yticks = np.linspace(f1,f2, num=7).astype(int)
dy_spec_yticks_y = np.linspace(0,ff-fs-1, len(dy_spec_yticks))

ax_2_0 = fig.add_subplot(g[2,0])
ax_2_0.imshow(P, cmap='Greys', aspect='auto', origin='lower', interpolation='none')
ax_2_0.set_xticks(dy_spec_xticks_x)
ax_2_0.set_xticklabels(dy_spec_xticks)
ax_2_0.set_yticks(dy_spec_yticks_y)
ax_2_0.set_yticklabels(dy_spec_yticks)
ax_2_0.set_xlabel('Time (ms)')
ax_2_0.set_ylabel('Frequency (MHz)')

# plot frequency autocorrelation
f_cor_xticks = np.round(np.linspace(-ms_tick, ms_tick, num=5),2)
f_cor_xticks_x = np.linspace(0,2*(nbin-1),num=len(f_cor_xticks))

ax_2_1 = fig.add_subplot(g[2,1])
ax_2_1.imshow(freq_corr, cmap='Reds', aspect='auto', origin='lower', interpolation='none')
ax_2_1.set_xticks(f_cor_xticks_x[1:-1])
ax_2_1.set_xticklabels(f_cor_xticks[1:-1])
ax_2_1.set_yticklabels([])
ax_2_1.set_yticks(dy_spec_yticks_y)	
ax_2_1.set_xlabel('Time Shift (ms)')

# plot phase autocorrelation
p_cor_yticks = np.linspace(-(f2-f1-1), f2-f1-1, num=7).astype(int)
p_cor_yticks_y = np.linspace(0,2*(nchan-1), len(p_cor_yticks))

ax_1_0 = fig.add_subplot(g[1,0])
ax_1_0.imshow(phase_corr, cmap='Blues', aspect='auto', origin='lower', interpolation='none')
ax_1_0.set_xticklabels([])
ax_1_0.set_xticks(dy_spec_xticks_x)
ax_1_0.set_yticks(p_cor_yticks_y[1:-1])
ax_1_0.set_yticklabels(p_cor_yticks[1:-1])
ax_1_0.set_ylabel('Frequency Shift (MHz)')

# plot 2D auto-correlation
ax_1_1 = fig.add_subplot(g[1,1])
ax_1_1.imshow(corr_2D, cmap='Purples', aspect='auto', origin='lower', interpolation='none')
ax_1_1.set_xticklabels([])
ax_1_1.set_xticks(f_cor_xticks_x[1:-1])
ax_1_1.set_yticklabels([])
ax_1_1.set_yticks(p_cor_yticks_y[1:-1])

## Summed auto-correlations
# plot summed frequency autocorrelation in frequency
ax_2_2 = fig.add_subplot(g[2,2])
ax_2_2.step(sum_freq_corr_freq, np.arange(len(sum_freq_corr_freq)), 
	color='red', where='mid', lw=1)
ax_2_2.set_ylim(0, nchan-1)
ax_2_2.set_xticklabels([])
ax_2_2.set_yticklabels([])
ax_2_2.set_yticks(dy_spec_yticks_y)
ax_2_2.set_xlabel('(arb. units)')

# plot summed frequency autocorrelation in time
ax_0_1 = fig.add_subplot(g[0,1])
ax_0_1.step(np.arange(len(sum_freq_corr_time)), sum_freq_corr_time, 
	color='red', where='mid', lw=1)
ax_0_1.set_xlim(0, 2*(nbin-1))
ax_0_1.set_xticklabels([])
ax_0_1.set_xticks(f_cor_xticks_x[1:-1])
ax_0_1.set_yticklabels([])

# plot summed phase autocorrelation in frequency
ax_1_2 = fig.add_subplot(g[1,2])
ax_1_2.step(sum_phase_corr_freq, np.arange(len(sum_phase_corr_freq)), 
	color='blue', where='mid', lw=1)
ax_1_2.set_ylim(0, 2*(nchan-1))
ax_1_2.set_xticklabels([])
ax_1_2.set_yticklabels([])
ax_1_2.set_yticks(p_cor_yticks_y[1:-1])

# plot summed phase autocorrelation in time
ax_0_0 = fig.add_subplot(g[0,0])
ax_0_0.step(np.arange(len(sum_phase_corr_time)), sum_phase_corr_time, 
	color='blue', where='mid', lw=1)
ax_0_0.set_xlim(0, nbin-1)
ax_0_0.set_xticklabels([])
ax_0_0.set_xticks(dy_spec_xticks_x)
ax_0_0.set_yticklabels([])
ax_0_0.set_ylabel('(arb. units)')

# plot summed 2D autocorrelation in frequency
ax_1_2.step(sum_corr_2D_freq, np.arange(len(sum_corr_2D_freq)), 
	color='purple', where='mid', lw=1)

# plot summed 2D autocorrelation in time
ax_0_1.step(np.arange(len(sum_corr_2D_time)), sum_corr_2D_time, 
	color='purple', where='mid', lw=1)

plt.savefig(sys.argv[0].split(os.extsep, 1)[0]+'_%s_'%sys.argv[2]+sys.argv[1].split(os.extsep, 1)[0]+'.pdf', dpi=300, bbox_inches='tight')
print(sys.argv[0].split(os.extsep, 1)[0]+'_%s_'%sys.argv[2]+sys.argv[1].split(os.extsep, 1)[0]+'.pdf')