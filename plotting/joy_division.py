import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

# Phase zoom factor
z = 0.005

# on-pulse phase start and finish
p1 = float(sys.argv[1])
p2 = float(sys.argv[2])

counter = 0
P = []
files = sorted(glob.glob("*.rescaled"))
for ar in files[0:10]:
	a = psrchive.Archive_load(ar)
	a.remove_baseline()
	a.tscrunch()
	a.fscrunch()
	a.pscrunch()
	data1 = a.get_data()
	nsub, npol, nchan, nbin = data1.shape



	# on-pulse phase bin start and finish
	ps = int(np.round(p1*nbin))
	pf = int(np.round(p2*nbin))

	profile = data1[0,0,0,ps:pf]

	P.append(profile)

	# file progress counter
	counter += 1
	print("%s/%s"%(counter,len(files))," files completed", end='\r')
# normalise
P = P/np.max(P)

#### PLOTTING ####
A4x, A4y = 8.27, 11.69
fig = plt.figure(figsize=(A4x, A4y),dpi=300)

# milliseconds per bin
period = a.integration_length()
bs = period/nbin
nbin_zoom = np.shape(P)[1]

### PLOT PULSE PROFILES
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))

m = 0.2
for i in range(np.shape(P)[0]):
	plt.plot(P[i]+i*m, linewidth=0.5, c='k')
	plt.fill_between(np.arange(pf-ps), P[i]+i*m, i*m, color='white')
	#plt.plot(np.arange(pf-ps), 0*np.arange(pf-ps)+i*200, linewidth=0.5, ls='--', color='r')

plt.xlim(0,pf-ps-1)
#plt.xticks(xticks_x, xticks)
plt.xlabel('Time (s)')
plt.ylabel('Pulse Number')

### SAVE FIGURE
plt.savefig('%s.pdf'%(sys.argv[0].split('.')[0]), bbox_inches='tight')
print('%s.pdf'%(sys.argv[0].split('.')[0]))
