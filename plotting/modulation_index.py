import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import os
import glob
from scipy.signal import peak_widths


def moments(n):
	"""Calculate moment of order n of spectra.
	"""
	m = []
	for i in range(nfile):
		m.append((1/nchan)*np.sum(S[i])**n)
	return m

def modulation_index():
	"""Calculate modulation index.
	"""
	MI = []
	# first and second moments
	m1 = moments(1)
	m2 = moments(2)
	for i in range(nfile):
		MI.append(np.sqrt((m2[i] - m1[i]**2)/m1[i]**2))
	return MI



# freq start and finish
f1 = float(sys.argv[1])
f2 = float(sys.argv[2])

# project code
PCODE = sys.argv[3]

# spectrum
if os.path.isfile("Pulse_Spectra_%s.npy"%PCODE):
	S = np.load("Pulse_Spectra_%s.npy"%PCODE)
else:
	S = []
	files = sorted(glob.glob("*.rescaled"))
	for ar in files:
		a = psrchive.Archive_load(ar)
		a.remove_baseline()
		a.tscrunch()
		a.pscrunch()
		data = a.get_data()
		nsub, npol, nchan, nbin = data.shape

		# peak and index
		pulse_profile = np.mean(data[0,0,:,:]/1000, axis=0)
		peak_idx = np.argmax(pulse_profile)

		# width of peaks for setting imshow widths
		width = peak_widths(pulse_profile, np.array([peak_idx]), rel_height=0.8)
		w = np.round(3*width[0]).astype(int)

		### FREQ ZOOM
		bw = a.get_bandwidth()
		cf = a.get_centre_frequency()
		#lowest observed frequency
		min_freq = cf-bw/2
		# fscrunching factor
		f_scr = bw/a.get_nchan()

		fs = int((f1-min_freq)/f_scr)
		ff = int((f2-min_freq)/f_scr)

		# spectrum
		# spectra start and finish
		s1 = np.round(width[2][0]).astype(int) 
		s2 = np.round(width[3][0]).astype(int)
		spectrum = np.mean(data[0,0,fs:ff,s1:s2]/1000, axis=1)
		S.append(spectrum)
	np.save("Pulse_Spectra_%s_%s_%s.npy"%(PCODE, int(f1), int(f2)), S)
nfile, nchan = np.shape(S)


### PLOTTING ###
if os.path.isfile("Modulation_Index_%s_%s_%s.npy"%(PCODE, int(f1), int(f2))):
	MI = np.load("Modulation_Index_%s_%s_%s.npy"%(PCODE, int(f1), int(f2)))
else:
	MI = modulation_index()
	np.save("Modulation_Index_%s_%s_%s.npy"%(PCODE, int(f1), int(f2)), MI)
print(MI)
A4x, A4y = 8.27, 11.69
fontsize = 10
fig = plt.figure(figsize=(A4x, A4x/2), dpi=300)
plt.hist(MI, bins=int(np.sqrt(len(MI))), edgecolor="k", color='white')
plt.axvline(np.mean(MI), color='r', linestyle='--')
plt.xlabel("Modulation Index", fontsize=fontsize)
plt.ylabel("Count", fontsize=fontsize)
plt.margins(x=0)
#plt.text(0.95, 0.95, 'MJD', transform=fig.transAxes, fontsize=fontsize, verticalalignment='top', horizontalalignment='right')
plt.savefig("Modulation_Index_%s_%s_%s.pdf"%(PCODE, int(f1), int(f2)), bbox_inches='tight', dpi=300)
print("Modulation_Index_%s_%s_%s.pdf"%(PCODE, int(f1), int(f2)))
