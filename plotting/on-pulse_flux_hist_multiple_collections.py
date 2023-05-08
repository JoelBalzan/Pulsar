import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob
import os

fluxes = np.empty(0)
counter = 0

if os.path.isfile('mean_fluxes.npy'):
	mean_fluxes = np.load('mean_fluxes.npy')
else:
	for ar in glob.glob("*.rescaled"):
		a = psrchive.Archive_load(ar)
		c1 = a.clone()
		c1.remove_baseline()
		c1.tscrunch()
		c1.fscrunch()
		c1.pscrunch()

		data1 = c1.get_data()
		nsub, npol, nchan, nbin = data1.shape

		# check paths for multiple collections
		path = os.readlink(ar)
		head = os.path.split(path)[0]

		if head == "/DATA/CARINA_6/bal256//J1809-1943/P970/37441/archive":
			bin_s = 0.55
			bin_f = 0.8
		elif head == "/DATA/CARINA_6/bal256//J1809-1943/PX500/38329/archive":
			bin_s = 0.27
			bin_f = 0.52
		elif head == "/DATA/CARINA_6/bal256//J1809-1943/PX500/38907/archive":
			bin_s = 0.68
			bin_f = 0.95
		elif head == "/DATA/CARINA_6/bal256//J1809-1943/PX500/39167/archive":
			bin_s = 0.1
			bin_f = 0.33

		# on-pulse start and finish phase bins
		on_s = int(bin_s*nbin)
		on_f = int(bin_f*nbin)

		# on-pulse mean flux density (Jy)
		flux = np.mean(data1[0,0,0,on_s:on_f]/1000, axis=0)
		mean_fluxes = np.concatenate((mean_fluxes, flux), axis=None)

		# file progress counter
		counter += 1
		print("%s/%s"%(counter,len(glob.glob("*.rescaled")))," files completed", end='\r')
	## SAVE FLUXES TO TXT FILE
	np.save('mean_fluxes.npy', mean_fluxes)


## MEAN OF ALL PULSES
avg = np.mean(mean_fluxes)
print("Mean = ", avg, " Jy")
print("Peak E = ", np.max(mean_fluxes)/avg, "<E>")
#fluxes_norm = mean_fluxes/avg


## PLOT HISTOGRAM
#hist, bins, _ = plt.hist(fluxes_norm, bins=20)
#logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

A4x, A4y = 8.27, 11.69
fontsize=20
plt.figure(figsize=(A4y,A4x),dpi=300)
#plt.rcParams["font.family"] = "serif"
#plt.rcParams["font.weight"] = "light"
ax = plt.subplot(111)
plt.hist(mean_fluxes, bins=20, edgecolor='black', color='white')
plt.axvline(avg, color='r', linewidth=1)
#ax.set_xscale("log")
#ax.set_yscale("log") 
plt.xlabel(r'$\langle{{E}}\rangle$ (Jy)', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
#plt.ylabel('log$_{10}$(Count)')
plt.ylabel('Count', fontsize=fontsize)
plt.yticks(fontsize=fontsize)
#plt.title('P970')

plt.savefig("%s_PALL.pdf"%sys.argv[0].split(os.extsep, 1)[0], bbox_inches='tight')
print("%s_PALL.pdf"%sys.argv[0].split(os.extsep, 1)[0])
