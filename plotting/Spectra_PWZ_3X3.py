import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import psrchive
import os
import glob
import string
from scipy.signal import find_peaks, peak_widths

# python *.py Polarisation freq_start freq_finish
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[1]
# freq start and finish
f1 = float(sys.argv[2])
f2 = float(sys.argv[3])

F = []
P = []
S = []
xticks = []
xticks_x = []
files = sorted(glob.glob("*.rescaled"))
# exit if <12 files
if len(files) < 12:
	print("12 files required. You only have %s files."%len(files))
	sys.exit()

for ar in files[0:len(files)]:
	if sys.argv[1] == "I":
		a = psrchive.Archive_load(ar)
		a.remove_baseline()
		a.tscrunch()
		a.pscrunch()
		data1 = a.get_data()
		nsub, npol, nchan, nbin = data1.shape

		# peak and index
		pulse_profile = np.mean(data1[0,0,:,:], axis=0)
		peak_idx = np.argmax(pulse_profile)

		# width of peaks for setting imshow widths
		w = np.round(3.4*peak_widths(pulse_profile, np.array([peak_idx]), rel_height=0.5)[0]).astype(int)

		# on-pulse phase bin start and finish
		ps = int(peak_idx - w)
		pf = int(peak_idx + w+1)
		# pulse profile
		F.append(pulse_profile[ps:pf]/1000)

		### FREQ ZOOM
		bw = a.get_bandwidth()
		cf = a.get_centre_frequency()
		#lowest observed frequency
		min_freq = cf-bw/2
		# fscrunching factor
		f_scr = bw/a.get_nchan()
	
		fs = int((f1-min_freq)/f_scr)
		ff = int((f2-min_freq)/f_scr)

		# intensity
		I = data1[0,0,fs:ff,ps:pf]/1000
		P.append(I)
		# spectrum
		spectrum = np.mean(I, axis=1)
		S.append(spectrum)

		# milliseconds per bin
		period = a.integration_length()
		bs = 1000*period/nbin
		nbin_zoom = np.shape(eval(p))[1]

		xt = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=5),2)#.astype(int)
		xticks.append(xt)
		xt_x = np.linspace(0,pf-ps-1,num=len(xt))
		xticks_x.append(xt_x)

	else:
		a = psrchive.Archive_load(ar)
		a.remove_baseline()
		a.tscrunch()
		a.bscrunch(8)
		data1 = a.get_data()
		nsub, npol, nchan, nbin = data1.shape

		# peak and index
		pulse_profile = data1.mean(axis=(1,2))[0]
		peak_idx = np.argmax(pulse_profile)

		# width of peaks for setting imshow widths
		w = np.round(3.4*peak_widths(pulse_profile, np.array([peak_idx]), rel_height=0.5)[0]).astype(int)

		# on-pulse phase bin start and finish
		ps = int(peak_idx - w)
		pf = int(peak_idx + w)
		# pulse profile
		F.append(pulse_profile[ps:pf]/1000)

		### FREQ ZOOM
		f_scr = (4032-704)/a.get_nchan()
		fs = int((f1-704)/f_scr)
		ff = int((f2-704)/f_scr)

		# polarisations
		if p == "SI":
			SI = data1[0,0,fs:ff,ps:pf]/1000
			P.append(SI)
			spectrum = np.mean(SI, axis=1)
			S.append(spectrum)
		if p == "SQ":
			SQ = data1[0,1,fs:ff,ps:pf]/1000
			P.append(SQ)
			spectrum = np.mean(SQ, axis=1)
			S.append(spectrum)
		if p == "SU":
			SU = data1[0,2,fs:ff,ps:pf]/1000
			P.append(SU)
			spectrum = np.mean(SU, axis=1)
			S.append(spectrum)
		if p == "L":
			SQ = data1[0,1,fs:ff,ps:pf]/1000
			SU = data1[0,2,fs:ff,ps:pf]/1000
			L = np.sqrt(SQ**2+SU**2)
			P.append(L)
			spectrum = np.mean(L, axis=1)
			S.append(spectrum)
		if p == "SV":
			SV = data1[0,3,fs:ff,ps:pf]/1000
			P.append(SV)
			spectrum = np.mean(SV, axis=1)
			S.append(spectrum)
		
		# milliseconds per bin
		period = a.integration_length()
		bs = 1000*period/nbin
		nbin_zoom = np.shape(eval(p))[1]

		xt = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=5),2)#.astype(int)
		xticks.append(xt)
		xt_x = np.linspace(0,pf-ps-1,num=len(xt))
		xticks_x.append(xt_x)


#### PLOTTING ####
if np.all(xticks == xticks[0]):
	bot = False
else:
	bot = True

A4x = 8.27
A4y = 11.69
fig = plt.figure(figsize=(A4x,A4y),dpi=600)
g = gridspec.GridSpec(ncols=3, nrows=4, hspace=0.11, wspace=0.11)

### PLOT DYNAMIC SPECTRA ###
yticks = np.linspace(f1,f2, num=7).astype(int)
yticks_y = np.linspace(0,ff-fs-1, len(yticks))

vmin = []
vmax = []
for i in range(len(files)):
	vmin.append(np.min(P[i]))
	vmax.append(0.8*np.max(P[i]))

# alphabet for plot labelling
alphabet = list(string.ascii_lowercase)

fontsize = 10
for i in range(len(files)):
	gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=g[i], hspace=0, wspace=0, height_ratios=[0.3,1], width_ratios=[1,0.3])
	ax00 = fig.add_subplot(gs[0,0])
	ax00.plot(F[i], color='k', linewidth=0.5)
	ax00.set_xticks(xticks_x[i])
	ax00.set_xticklabels([])
	ax00.set_yticks([])
	ax00.set_yticklabels([])
	ax00.margins(x=0)

	ax11 = fig.add_subplot(gs[1,1])
	ax11.plot(S[i], np.arange(ff-fs), color='k', linewidth=0.5)
	ax11.set_xticks([])
	ax11.set_xticklabels([])
	ax11.set_yticks(yticks_y)
	ax11.set_yticklabels([])
	ax11.set_ylim(0,ff-fs-1)

	ax10 = fig.add_subplot(gs[1,0])
	ax10.imshow(P[i], cmap="viridis", 
			  vmin=vmin[i], 
			  vmax=vmax[i], 
			  aspect='auto', origin='lower', interpolation='none')
	ax10.set_xticks(xticks_x[i])
	ax10.set_xticklabels(xticks[i], fontsize=fontsize)
	plt.yticks(yticks_y, yticks, fontsize=fontsize)
	if i == 0 or i == 3 or i == 6 or i == 9:
		ax10.set_ylabel('Frequency (MHz)', fontsize=fontsize)
		ax10.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=True, right=True, top=True)
		ax10.text(0.05, 0.95, alphabet[i]+")", transform=ax10.transAxes, fontsize=fontsize, fontweight='bold', va='top', color='w')
	if i == 1 or i == 2 or i == 4 or i == 5 or i == 7 or i == 8:
		ax10.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=False, right=True, top=True)
		ax10.text(0.05, 0.95, alphabet[i]+")", transform=ax10.transAxes, fontsize=fontsize, fontweight='bold', va='top', color='w')
	if i == 9:
		ax10.set_xlabel('Time (s)', fontsize=fontsize)
		ax10.set_ylabel('Frequency (MHz)', fontsize=fontsize)
		ax10.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True, right=True, top=True)
	if i == 10 or i == 11:
		ax10.set_xlabel('Time (s)', fontsize=fontsize)
		ax10.tick_params(bottom=True, labelbottom=True, left=True, labelleft=False, right=True, top=True)
		ax10.text(0.05, 0.95, alphabet[i]+")", transform=ax10.transAxes, fontsize=fontsize, fontweight='bold', va='top', color='w')


### SAVE FIGURE
plt.savefig('Spectra_PWZ_%s_PX500.png'%p, bbox_inches='tight')
print('Spectra_PWZ_%s_PX500.png'%p)
