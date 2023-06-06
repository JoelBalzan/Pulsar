import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import psrchive



def subband_flux(DS, nsub, on_s, on_f):
	for i in range(nsub):
	    for j in range(n_subbands):
	    	subbands["sband"+str(j)].append(np.sum(np.mean(DS[i,0,j*n_chan_profile:(j+1)*n_chan_profile,on_s:on_f],axis=0))/(on_f-on_s))

# import folded files
P970 = psrchive.Archive_load('P970_fold.rf')
P970.pscrunch()
P970_data = P970.get_data()
nsub_P970, npol, nchan, nbin = P970_data.shape
P970_mjds = P970.get_mjds()
on_s_P970 = int(0.53*nbin)
on_f_P970 = int(0.78*nbin)

PX500_38329 = psrchive.Archive_load('PX500_38329_fold.rf')
PX500_38329.pscrunch()
PX500_38329_data = PX500_38329.get_data()
nsub_PX500_38329, npol, nchan, nbin = PX500_38329_data.shape
PX500_38329_mjds = PX500_38329.get_mjds()
on_s_PX500_38329 = int(0.27*nbin)
on_f_PX500_38329 = int(0.51*nbin)

PX500_38907 = psrchive.Archive_load('PX500_38907_fold.rf')
PX500_38907.pscrunch()
PX500_38907_data = PX500_38907.get_data()
nsub_PX500_38907, npol, nchan, nbin = PX500_38907_data.shape
PX500_38907_mjds = PX500_38907.get_mjds()
on_s_PX500_38907 = int(0.68*nbin)
on_f_PX500_38907 = int(0.94*nbin)

#PX500_39167 = psrchive.Archive_load('PX500_39167_fold.rf')
#PX500_39167.pscrunch()
#PX500_39167_data = PX500_39167.get_data()
#nsub_PX500_39167, npol, nchan, nbin = PX500_39167_data.shape
#PX500_39167_mjds = PX500_39167.get_mjds()
#on_s_PX500_39167 = int(0.11*nbin)
#on_f_PX500_39167 = int(0.33*nbin)

n_subbands = 13
n_chan_profile = int(nchan/n_subbands)

subbands={}
for i in range(n_subbands):
	key = str("sband"+str(i))
	subbands[key] = []
#for key,value in subbands.items():
#	exec(f'{key}={value}')
subband_flux(P970_data, nsub_P970, on_s_P970, on_f_P970)
subband_flux(PX500_38329_data, nsub_PX500_38329, on_s_PX500_38329, on_f_PX500_38329)
subband_flux(PX500_38907_data, nsub_PX500_38907, on_s_PX500_38907, on_f_PX500_38907)
#subband_flux(PX500_39167_data, nsub_PX500_39167, on_s_PX500_39167, on_f_PX500_39167)
np.save('subband_flux.npy',subbands, allow_pickle=True)

A4x, A4y = 8.27, 11.69
fig = plt.figure(figsize=(A4y, A4x), dpi=300)

for i in range(n_subbands):
	plt.plot(subbands["sband"+str(i)],label='sband'+str(i))
plt.legend()
plt.savefig(sys.argv[0].split(os.extsep, 1)[0]+'.png', dpi=300, bbox_inches='tight')
