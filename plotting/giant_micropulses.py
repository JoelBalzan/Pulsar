import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob
import os


if os.path.isfile('peak_flux.npy'):
	peak_flux = np.load('peak_flux.npy')
else:
	peak_flux = []
	counter = 0
	for ar in glob.glob("*.rescaled"):
		a = psrchive.Archive_load(ar)
		c1 = a.clone()
		c1.remove_baseline()
		c1.tscrunch()
		c1.fscrunch()
		c1.pscrunch()

		data1 = c1.get_data()
		nsub, npol, nchan, nbin = data1.shape

		# on-pulse start and finish phase bins
		on_s = int(float(sys.argv[1])*nbin)
		on_f = int(float(sys.argv[2])*nbin)

		# on-pulse mean flux density (Jy)
		flux = np.mean(data1[0,0,0,on_s:on_f]/1000, axis=0)
		peak_flux.append(flux)

		# file progress counter
		counter += 1
		print("%s/%s"%(counter,len(glob.glob("*.rescaled")))," files completed", end='\r')
	## SAVE FLUXES TO TXT FILE
	peak_flux = np.array(peak_flux)
	np.save('fluxes.txt', peak_flux)

## MEAN OF ALL PULSES
avg = np.mean(peak_flux)
print("Mean = ", avg, " Jy")
fluxes_norm = peak_flux/avg


## PLOT HISTOGRAM
#hist, bins, _ = plt.hist(fluxes, bins=20)
#logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

plt.figure(figsize=(15,10),dpi=300)
ax = plt.subplot(111)
plt.hist(peak_flux, bins=20, edgecolor='black', color='white')
plt.axvline(avg, color='r', linewidth=1)
#ax.set_xscale("log")
#ax.set_yscale("log") 
plt.xlabel('Peak Flux Density (Jy)')
#plt.ylabel('log$_{10}$(Count)')
plt.ylabel('Count')
#plt.title('P970')

plt.savefig("%s_PALL.pdf"%sys.argv[0].split(os.extsep, 1)[0])
print("%s_PALL.pdf"%sys.argv[0].split(os.extsep, 1)[0])
