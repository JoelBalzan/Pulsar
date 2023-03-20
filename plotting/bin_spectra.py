import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.fft
from matplotlib import gridspec
import psrchive
import os
from scipy.signal import find_peaks, peak_widths

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

a = psrchive.Archive_load(sys.argv[1])
# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[2]


### FREQ ZOOM
f_scr = (4032-704)/a.get_nchan()

f1 = float(sys.argv[3])
f2 = float(sys.argv[4])
fs = int((f1-704)/f_scr)
ff = int((f2-704)/f_scr)


if sys.argv[2] == "I":
	a.remove_baseline()
	a.tscrunch()
	a.pscrunch()
	data2 = a.get_data()
	nsub, npol, nchan, nbin = data2.shape

    # peak and index
	# flux in Jy
	flux = np.mean(data2[0,0,:,:], axis=0)/1000
	peak_idx = np.array([np.argmax(flux)])
	w = np.round(peak_widths(flux, peak_idx, rel_height=0.5)).astype(int)
	
	# intensity
	I = data2[0,0,fs:ff,w[2][0]:w[3][0]]
	w[2] -= 1
	w[3] += 1
	
else:
	a.remove_baseline()
	a.tscrunch()
	#a.bscrunch(8)
	data = a.get_data()
	nsub, npol, nchan, nbin = data.shape

	# peak and index
	# flux in Jy
	flux = data.mean(axis=(1,2))[0]
	peak_idx = np.array([np.argmax(flux)])
	w = np.round(peak_widths(flux, peak_idx, rel_height=0.8)).astype(int)
	w[2] -= 1
	w[3] += 1

	# polarisations
	if p == "SI":
		SI = data[0,0,fs:ff,w[2][0]:w[3][0]]
	if p == "SQ":
		SQ = data[0,1,fs:ff,w[2][0]:w[3][0]]
	if p == "SU":
		SU = data[0,2,fs:ff,w[2][0]:w[3][0]]
	if p == "L":
		L = np.sqrt(data[0,1,fs:ff,w[2][0]:w[3][0]]**2+data[0,2,fs:ff,w[2][0]:w[3][0]]**2)
	if p == "SV":
		SV = data[0,3,fs:ff,w[2][0]:w[3][0]]

# seconds per bin
nbin_zoom = np.shape(eval(p))[1]

### DEFINE SPECTRA 
S = []
for i in range(nbin_zoom):
	S.append(eval(p)[:,i])
S = np.array(S)
# normalise spectra
S = S/S.max()
#S[S == 0.] = np.nan

#### PLOTTING ####
fig = plt.figure(figsize=(15,10),dpi=300)

### PLOT SPECTRA
xticks = np.linspace(f1,f2, num=14).astype(int)
xticks_x = np.linspace(0,ff-fs-1, len(xticks))

for i in range(nbin_zoom):
	plt.plot(S[i]-i/2, linewidth=0.5, c='k')
	plt.plot(np.arange(ff-fs), 0*np.arange(ff-fs)-i/2, linewidth=0.5, ls='--', color='r')
	
plt.xticks(xticks_x, xticks)
plt.yticks([])
plt.xlabel('Frequency (MHz)')
plt.ylabel('Phase Bins')

#plt.title('%s Polarisation %s'%(p,sys.argv[1].split('.')[0]))



### SAVE FIGURE
plt.savefig('bin_spectra_%s_%s_%s_%s.pdf'%(p, sys.argv[1].split(os.extsep, 1)[0], int(f1), int(f2)), bbox_inches='tight')
print('bin_spectra_%s_%s_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0], int(f1), int(f2)))



### PLOT FFT
#fig = plt.figure(figsize=(15,10),dpi=300)
#
## manually select bins from previous plot
##S[S == np.nan] = 0.
#
#FT = []
#for i in (np.arange(7,10,1)):
#	FT.append(scipy.fft.fft(S[i]).real)
#FT = np.array(FT)
#FT = FT/FT.max()
#
#for i in range(np.shape(FT)[0]):
#	plt.plot(FT[i,1:np.shape(FT)[1]//2]-i/2, linewidth=0.5, c='k')
#	plt.plot(np.arange(ff-fs), 0*np.arange(ff-fs)-i/2, linewidth=0.5, ls='--', color='r')
#
#plt.xlabel('Flux')
#plt.ylabel('Normalised Intensity')
#
#### SAVE FIGURE
#plt.savefig('fft_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
#print('fft_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))