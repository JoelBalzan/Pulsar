import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.fft
from matplotlib import gridspec
import psrchive
import os
from scipy.signal import find_peaks

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

a = psrchive.Archive_load(sys.argv[1])
# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[2]


### DETERMINE PEAK FLUX AND INDEX FOR PLOT CENTRING
c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape

# neutron star period (ms)
period = a.integration_length()


# peak and index
flux = data[0,0,0,:]
peaks, _ = find_peaks(flux)
peak_flux = np.sort(flux[peaks])[-1]
peak_idx = np.where(flux==peak_flux)[0][0]

# on-pulse phase start and finish
#p1 = np.round(peak_idx/nbin - 0.1, 4)
#p2 = np.round(peak_idx/nbin + 0.1, 4)
z = 0.0002
p1 = np.round(peak_idx/nbin - z, 4)
p2 = np.round(peak_idx/nbin + z, 4)


### FREQ ZOOM
f_scr = (4032-704)/a.get_nchan()

f1 = float(sys.argv[3])
f2 = float(sys.argv[4])
fs = int((f1-704)/f_scr)
ff = int((f2-704)/f_scr - 1) #extra -1 to account for indexing


if sys.argv[2] == "I":
	c1 = a.clone()
	c1.remove_baseline()
	c1.tscrunch()
	c1.pscrunch()
	data2 = c1.get_data()
	nsub, npol, nchan, nbin = data2.shape

	# on-pulse phase bin start and finish
	ps = int(np.round(p1*nbin))
	pf = int(np.round(p2*nbin))

	# intensity
	I = data2[0,0,fs:ff,ps:pf]
    
else:
	c1 = a.clone()
	c1.remove_baseline()
	c1.tscrunch()
	c1.bscrunch(8)
	data2 = c1.get_data()
	nsub, npol, nchan, nbin = data2.shape

	# on-pulse phase bin start and finish
	ps = int(np.round(p1*nbin))
	pf = int(np.round(p2*nbin))

	# polarisations
	SI = data2[0,0,fs:ff,ps:pf]
	SQ = data2[0,1,fs:ff,ps:pf]
	SU = data2[0,2,fs:ff,ps:pf]
	L = np.sqrt(data2[0,1,fs:ff,ps:pf]**2+data2[0,2,fs:ff,ps:pf]**2)
	SV = data2[0,3,fs:ff,ps:pf]

# seconds per bin
bs = 1000*period/nbin
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
plt.xlabel('Frequency (MHz)')
plt.ylabel('Bin')

plt.title('%s Polarisation %s'%(p,sys.argv[1].split('.')[0]))



### SAVE FIGURE
plt.savefig('bin_spectra_%s_%s_%s_%s.pdf'%(p, sys.argv[1].split(os.extsep, 1)[0], int(f1), int(f2)))
print('bin_spectra_%s_%s_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0], int(f1), int(f2)))



### PLOT FFT
fig = plt.figure(figsize=(15,10),dpi=300)

# manually select bins from previous plot
#S[S == np.nan] = 0.

FT = []
for i in (np.arange(7,10,1)):
	FT.append(scipy.fft.fft(S[i]).real)
FT = np.array(FT)
FT = FT/FT.max()

for i in range(np.shape(FT)[0]):
	plt.plot(FT[i,1:np.shape(FT)[1]//2]-i/2, linewidth=0.5, c='k')
	plt.plot(np.arange(ff-fs), 0*np.arange(ff-fs)-i/2, linewidth=0.5, ls='--', color='r')

plt.xlabel('Flux')
plt.ylabel('Normalised Intensity')

### SAVE FIGURE
plt.savefig('fft_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
print('fft_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))