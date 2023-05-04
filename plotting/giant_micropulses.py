import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob
import os

PCODE = "PX500_39167"
if os.path.isfile('peak_flux_'+PCODE+'.npy'):
	peak_flux = np.load('peak_flux_'+PCODE+'.npy')
else:
	peak_flux = []
	counter = 0
	for ar in glob.glob("*.rescaled"):
		a = psrchive.Archive_load(ar)
		a.remove_baseline()
		a.tscrunch()
		a.fscrunch()
		a.pscrunch()
		data = a.get_data()
		peak_flux.append(np.max(data[0,0,0,:])/1000)

		# file progress counter
		counter += 1
		print("%s/%s"%(counter,len(glob.glob("*.rescaled")))," files completed", end='\r')
	## SAVE FLUXES TO FILE
	peak_flux = np.array(peak_flux)
	np.save('peak_flux_'+PCODE+'.txt', peak_flux)

## MEAN OF ALL PULSES
avg = np.mean(peak_flux)
print("Mean = ", avg, " Jy")
print("Maximum =", np.max(peak_flux), " Jy =", np.max(peak_flux)/avg, " times the mean")
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

plt.savefig("%s_%s.pdf"%(sys.argv[0].split(os.extsep, 1)[0], PCODE), bbox_inches='tight')
print("%s_%s.pdf"%(sys.argv[0].split(os.extsep, 1)[0], PCODE))
