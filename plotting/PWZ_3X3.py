import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import psrchive
import os
import glob
import string
from scipy.signal import peak_widths, find_peaks

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[1]
# freq start and finish
f1 = float(sys.argv[2])
f2 = float(sys.argv[3])

# dynamic spectrum
P = []
xticks = []
xticks_x = []
files = sorted(glob.glob("*.rescaled"))
# exit if <12 files
if len(files) < 12:
	print("12 files required. You only have %s files."%len(files))
	sys.exit()

for ar in files:
	if sys.argv[1] == "I":
		a = psrchive.Archive_load(ar)
		a.remove_baseline()
		a.tscrunch()
		a.pscrunch()
		data1 = a.get_data()
		nsub, npol, nchan, nbin = data1.shape

		# peak and index
		profile = np.mean(data1[0,0,:,:], axis=0)
		peak_idx = np.argmax(profile)
		if ar == 'pulse_65080037.calib.rescaled':
			peaks, _ = find_peaks(profile)
			peak_idx = np.where(profile==np.sort(profile[peaks])[-2])[0][0]
		
		# width of peaks for setting imshow widths
		width = peak_widths(profile, np.array([peak_idx]), rel_height=0.8)
		w = np.round(3*width[0]).astype(int)

		# on-pulse phase bin start and finish
		#ps = int(peak_idx - w)
		#pf = int(peak_idx + w + 1)
		ps = int(peak_idx - 12)
		pf = int(peak_idx + 12 + 1)

		### FREQ ZOOM
		bw = a.get_bandwidth()
		cf = a.get_centre_frequency()
		#lowest observed frequency
		min_freq = cf-bw/2
		# fscrunching factor
		f_scr = (4032-704)/a.get_nchan()

		fs = int((f1-704)/f_scr)
		ff = int((f2-704)/f_scr)

		# intensity
		I = data1[0,0,fs:ff,ps:pf]

		P.append(I)

		# milliseconds per bin
		period = a.integration_length()
		bs = 1000*period/nbin
		nbin_zoom = np.shape(eval(p))[1]

		xt = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=5),1)[1:-1]#.astype(int)
		xticks.append(xt)
		xt_x = np.linspace(0,pf-ps-1,num=len(xt)+2)[1:-1]
		xticks_x.append(xt_x)

	else:
		a = psrchive.Archive_load(ar)
		a.remove_baseline()
		a.tscrunch()
		a.bscrunch(8)
		data1 = a.get_data()
		nsub, npol, nchan, nbin = data1.shape

		# peak and index
		profile = data1.mean(axis=(1,2))[0]
		peak_idx = np.argmax(profile)
		width = peak_widths(profile, np.array([peak_idx]), rel_height=0.8)
		# width in bins
		w = np.round(3.4*width[0]).astype(int)

		# on-pulse phase bin start and finish
		ps = int(peak_idx - w)
		pf = int(peak_idx + w+1)

		### FREQ ZOOM
		f_scr = (4032-704)/a.get_nchan()

		fs = int((f1-704)/f_scr)
		ff = int((f2-704)/f_scr)

		# polarisations
		if p == "SI":
			SI = data1[0,0,fs:ff,ps:pf]
			P.append(SI)
		if p == "SQ":
			SQ = data1[0,1,fs:ff,ps:pf]
			P.append(SQ)
		if p == "SU":
			SU = data1[0,2,fs:ff,ps:pf]
			P.append(SU)
		if p == "L":
			L = np.sqrt(data1[0,1,fs:ff,ps:pf]**2+data1[0,2,fs:ff,ps:pf]**2)
			P.append(L)
		if p == "SV":
			SV = data1[0,3,fs:ff,ps:pf]
			P.append(SV)

		# milliseconds per bin
		period = a.integration_length()
		bs = 1000*period/nbin
		nbin_zoom = np.shape(eval(p))[1]

		xt = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=5),2)[1:-1]#.astype(int)
		xticks.append(xt)
		xt_x = np.linspace(0,pf-ps-1,num=len(xt)+2)[1:-1]
		xticks_x.append(xt_x)

#### PLOTTING ####
if np.all(xticks == xticks[0]):
	bot = False
else:
	bot = True

A4x = 8.27
A4y = 11.69
fig = plt.figure(figsize=(A4x,A4y/1.5),dpi=600)
g = gridspec.GridSpec(ncols=3, nrows=4, hspace=0, wspace=0)


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
	ax = fig.add_subplot(g[i])
	ax.imshow(P[i], cmap="Greys", 
			  vmin=vmin[i], 
			  vmax=vmax[i], 
			  aspect='auto', origin='lower', interpolation='none')
	ax.set_xticks(xticks_x[i])
	ax.set_xticklabels(xticks[i], fontsize=fontsize)
	plt.yticks(yticks_y, yticks, fontsize=fontsize)
	if i == 0 or i == 3 or i == 6 or i == 9:
		ax.set_ylabel('Frequency (MHz)', fontsize=fontsize)
		ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=True, right=True, top=True)
		ax.text(0.05, 0.95, alphabet[i]+")", transform=ax.transAxes, fontsize=fontsize, fontweight='bold', va='top', color='w')
	if i == 1 or i == 2 or i == 4 or i == 5 or i == 7 or i == 8:
		ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=False, right=True, top=True)
		ax.text(0.05, 0.95, alphabet[i]+")", transform=ax.transAxes, fontsize=fontsize, fontweight='bold', va='top', color='w')
	if i == 9:
		ax.set_xlabel('Time (s)', fontsize=fontsize)
		ax.set_ylabel('Frequency (MHz)', fontsize=fontsize)
		ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True, right=True, top=True)
	if i == 10 or i == 11:
		ax.set_xlabel('Time (s)', fontsize=fontsize)
		ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=False, right=True, top=True)
		ax.text(0.05, 0.95, alphabet[i]+")", transform=ax.transAxes, fontsize=fontsize, fontweight='bold', va='top', color='w')
	if i == 0 or i == 6:
		ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=False, right=True, top=True)
		ax.set_ylabel("")



### SAVE FIGURE
plt.savefig('PWZ_%s_PX500.png'%p, bbox_inches='tight')
print('PWZ_%s_PX500.png'%p)
