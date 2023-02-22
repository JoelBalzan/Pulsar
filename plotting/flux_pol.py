import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import psrchive


a = psrchive.Archive_load(sys.argv[1])
c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape


L = np.sqrt(np.power(data[0,1,0,:], 2.) + np.power(data[0,2,0,:], 2.))
V = data[0,3,0,:]


p1 = float(sys.argv[2])
p2 = float(sys.argv[3])
ps = int(p1*nbin)
pf = int(p2*nbin)

ticks = np.round(np.linspace(p1,p2,num=11),2)
ticks_x = np.linspace(0,pf-ps-1,num=len(ticks))

plt.figure(figsize=(15,10),dpi=300)
ax = plt.subplot(111)
ax.set_xlim(0,pf-ps)
plt.plot(L[ps:pf], zorder=2, label='L', c='r')
plt.plot(V[ps:pf], zorder=3, label='V', c='b')


c1 = a.clone()
c1.remove_baseline()
c1.tscrunch()
c1.fscrunch()
c1.pscrunch()
data1 = c1.get_data()
nsub, npol, nchan, nbin = data1.shape


plt.plot(data1[0,0,0,ps:pf], zorder=1, label='Flux Density', c='black')
plt.xticks(ticks_x, ticks)
plt.xlabel('Phase')
plt.ylabel('Flux Density (mJy)')
plt.title('%s'%sys.argv[1].split('.')[0])
plt.legend()

plt.savefig("flux_pol_%s.pdf"%sys.argv[1].split('.')[0])
print("flux_pol_%s.pdf"%sys.argv[1].split('.')[0])