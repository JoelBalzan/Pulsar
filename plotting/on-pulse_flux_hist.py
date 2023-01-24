import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob
import os

fluxes = np.empty(0)
counter = 0

if os.path.isfile('fluxes.txt'):
	fluxes = np.loadtxt('fluxes.txt')
else:
	for ar in glob.glob("*.rescaled"):
		a = psrchive.Archive_load(ar)
		#c = a.clone()
		#c.tscrunch()
		#data = c.get_data()
		#nsub, npol, nchan, nbin = data.shape
	#
	#
		#delta = float(sys.argv[1])
		#I = data[0,0,:,:]
		#I = np.nanmean(I, axis=0)
		#num = int(nbin*delta)
		#m = nbin
	#
		#std = []
		#for i in range(m):
		#	if (i+num) < m:
		#		temp = I[i:(i+num)]
		#		std.append(np.nanstd(temp))
		#	else:
		#		temp = np.concatenate((I[i:m],I[0:(i+num-m)]))
		#		std.append(np.nanstd(temp))
	#
		#std = np.array(std)
		#on_s = int(np.argmax(std))
		#on_f = int(on_s + 0.25*nbin)
	#
		# load data
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
		flux = np.sum(data1[0,0,0,on_s:on_f]/1000)/(on_f-on_s)
		fluxes = np.concatenate((fluxes, flux), axis=None)

		# file progress counter
		counter = counter + 1
		print("%s/%s"%(counter,len(glob.glob("*.rescaled")))," files completed", end='\r')
	## SAVE FLUXES TO TXT FILE
	np.savetxt('fluxes.txt', fluxes)

## MEAN OF ALL PULSES
avg = np.average(fluxes)
print("Mean = ", avg, " Jy")
fluxes_norm = fluxes/avg


## PLOT HISTOGRAM
hist, bins, _ = plt.hist(fluxes_norm, bins=20)
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

plt.figure(figsize=(15,10),dpi=300)
ax = plt.subplot(111)
plt.hist(fluxes_norm, bins=logbins, edgecolor='black', color='white')
#plt.axvline(avg, color='k', linestyle='dashed', linewidth=1)
ax.set_xscale("log")
#ax.set_yscale("log") 
plt.xlabel('log$_{10}$(Flux/Mean Flux)')
#plt.ylabel('log$_{10}$(Count)')
plt.ylabel('Count')
#plt.title('P970')

plt.savefig("%s_PALL_log_norm.pdf"%sys.argv[0].split(os.extsep, 1)[0])
print("%s_PALL_log_norm.pdf"%sys.argv[0].split(os.extsep, 1)[0])
