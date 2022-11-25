import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive


a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.fscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape


L = np.sqrt(np.power(data[0,1,0,:], 2.) + np.power(data[0,2,0,:], 2.))
V = data[0,3,0,:]

x = np.arange(nbin)

p1 = float(sys.argv[2])
p2 = float(sys.argv[3])
ps = int(p1*nbin)
pf = int(p2*nbin)

x = np.arange(pf-ps)
ticks = np.round(np.linspace(p1,p2,num=11),2)
ticks_x = np.linspace(0,pf-ps,num=11)

plt.figure(figsize=(15,10),dpi=300)
plt.plot(x, L[ps:pf], zorder=2, label='L', c='r')
plt.plot(x, V[ps:pf], zorder=3, label='V', c='b')

a.pscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

plt.plot(x, data[0,0,0,ps:pf], zorder=1, label='Flux Density', c='black')
plt.xticks(ticks_x, ticks)
plt.xlabel('Phase')
plt.ylabel('Flux Density')
plt.legend()

plt.savefig("FP.png")
