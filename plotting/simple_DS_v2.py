import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import os
import glob
from scipy.signal import find_peaks, peak_widths

# python *.py PCODE freq_start freq_finish phase_start phase_finish
# project code
PCODE = sys.argv[1]

# freq start and finish
f1 = float(sys.argv[2])
f2 = float(sys.argv[3])

files = sorted(glob.glob("pulse_*.calib.rescaled"))
# period in minutes
f = psrchive.Archive_load(files[0])
period = f.integration_length()/60

### FREQ ZOOM
bw = f.get_bandwidth()
cf = f.get_centre_frequency()
#lowest observed frequency
min_freq = cf-bw/2
# fscrunching factor
f_scr = bw/f.get_nchan()
fs = int((f1-min_freq)/f_scr)
ff = int((f2-min_freq)/f_scr)

### PHASE ZOOM
nbins = f.get_nbin()
del f

# minimum height of peaks (Jy)
h=2

# spectrum
if os.path.isfile("%s_Spectra_v2_704_4032_on_%s.npy"%(PCODE, h)):
	S = np.load("%s_Spectra_v2_704_4032_on_%s.npy"%(PCODE, h))
	S = S[fs:ff, :]
else:
	p1 = int(float(sys.argv[4])*nbins)
	p2 = int(float(sys.argv[5])*nbins)
	S = []
	for ar in files:
		a = psrchive.Archive_load(ar)
		a.remove_baseline()
		a.tscrunch()
		a.pscrunch()
		data = a.get_data()
		data = np.squeeze(data[0,0,:,p1:p2])
		nchan, nbin = data.shape

		# pulse profile in Jy
		profile = np.mean(data[:,:], axis=0)/1000

		peaks, _ = find_peaks(profile, height=h, distance=8)
		#peak widths
		width = peak_widths(profile, peaks, rel_height=0.8)

		peak_spectra = []
		for i in range(len(peaks)):
			s1 = np.round(width[2][i]).astype(int) - 1
			s2 = np.round(width[3][i]).astype(int) + 1
			s = np.mean(data[:,s1:s2], axis=1)
			peak_spectra.append(s)
		S.append(np.mean(peak_spectra, axis=0))
	S = np.rot90(S, k=3)
	np.save("%s_Spectra_v2_%s_%s_on_%s.npy"%(PCODE, int(f1), int(f2), h), S)
nfile, nchan = np.shape(S)


### PLOTTING ###
A4x, A4y = 8.27, 11.69
fontsize = 10



# axes ticks
m_tick = nfile*period
xticks = np.round(np.linspace(0,m_tick,num=9),2)
xticks_x = np.linspace(0,nfile-1,num=len(xticks))
yticks = np.linspace(f1,f2, num=13).astype(int)
yticks_y = np.linspace(0,ff-fs-1, len(yticks))

# Choose good min/max values to have a good dynamic range
mn = np.mean(S)
std = np.std(S)
vmin = mn - 3*std
vmax = mn + 7*std

fig = plt.figure(figsize=(A4x, A4x), dpi=300)
plt.imshow(S, cmap='inferno', aspect='auto', interpolation='none', origin='lower', vmin=vmin, vmax=vmax)
plt.xticks(xticks_x, xticks, fontsize=fontsize)
plt.xlabel("Time (m)", fontsize=fontsize)
plt.yticks(yticks_y, yticks, fontsize=fontsize)
plt.ylabel("Frequency (MHz)", fontsize=fontsize)
plt.margins(x=0)

plt.savefig("simple_DS_v2_%s_%s_%s_on_%s.pdf"%(PCODE, int(f1), int(f2), h), bbox_inches='tight', dpi=300)
print("simple_DS_v2_%s_%s_%s_on_%s.pdf"%(PCODE, int(f1), int(f2), h))