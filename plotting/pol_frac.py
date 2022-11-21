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


I = data[0,0,0,:]
I[I==0] = np.nan
Q = data[0,1,0,:]
U = data[0,2,0,:]
V = data[0,3,0,:]

l = np.float64(Q+U)/I
c = np.float64(V)/I

p1 = float(sys.argv[2])
p2 = float(sys.argv[3])
ps = int(p1*nbin)
pf = int(p2*nbin)

x = np.arange(pf-ps)
ticks = np.round(np.arange(p1,p2+0.1, step=0.1), 1)
ticks_x = np.arange(0,pf-ps, step=int((pf-ps)/(len(ticks)-1)))

plt.figure(figsize=(15,10),dpi=300)
plt.plot(x, l[ps:pf], label='L/I')
plt.plot(x, c[ps:pf], label='V/I')
plt.xticks(ticks_x, ticks)
plt.xlabel('Phase')
plt.ylabel('Polarisation Fraction')
plt.legend()

plt.savefig("pf.png")