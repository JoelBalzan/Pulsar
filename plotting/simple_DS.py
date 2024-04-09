import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import os
import glob


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

# spectrum
if os.path.isfile("%s_Spectra_704_4032_on.npy"%(PCODE)):
	S = np.load("%s_Spectra_704_4032_on.npy"%(PCODE))
	S = S[fs:ff, :]
else:
	S = []
	p1 = float(sys.argv[4])*nbins
	p2 = float(sys.argv[5])*nbins
	for ar in files:
		a = psrchive.Archive_load(ar)
		a.remove_baseline()
		a.tscrunch()
		a.pscrunch()
		data = a.get_data()
		nsub, npol, nchan, nbin = data.shape

		spectrum = np.mean(data[0,0,fs:ff,int(p1):int(p2)], axis=1)
		del a

		S.append(spectrum)
		del spectrum
	S = np.rot90(S, k=3)
	np.save("%s_Spectra_%s_%s_on.npy"%(PCODE, int(f1), int(f2)), S)
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

fig = plt.figure(figsize=(A4x, A4x), dpi=300)
plt.imshow(S, cmap='inferno', aspect='auto', interpolation='none', origin='lower', vmin=0.1*np.min(S), vmax=0.7*np.max(S))
plt.xticks(xticks_x, xticks, fontsize=fontsize)
plt.xlabel("Time (m)", fontsize=fontsize)
plt.yticks(yticks_y, yticks, fontsize=fontsize)
plt.ylabel("Frequency (MHz)", fontsize=fontsize)
plt.margins(x=0)

plt.savefig("simple_DS_%s_%s_%s_on.pdf"%(PCODE, int(f1), int(f2)), bbox_inches='tight', dpi=300)
print("simple_DS_%s_%s_%s_on.pdf"%(PCODE, int(f1), int(f2)))

