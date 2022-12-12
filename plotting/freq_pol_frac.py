import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import pandas as pd
import os

a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.fscrunch(4)
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

# on pulse phase start
p1 = float(sys.argv[2]) 
# on pulse phase finish
p2 = float(sys.argv[3])
ps = int(p1*nbin)
pf = int(p2*nbin)


I = []
for i in range(nchan):
    s1 = pd.DataFrame([[np.average(data[0,0,i,ps:pf])]])
    I.append(s1)
I = pd.concat(I, ignore_index=True)

Q = []
for i in range(nchan):
    s2 = pd.DataFrame([[np.average(data[0,1,i,ps:pf])]])
    Q.append(s2)
Q = pd.concat(Q, ignore_index=True)

U = []
for i in range(nchan):
    s3 = pd.DataFrame([[np.average(data[0,2,i,ps:pf])]])
    U.append(s3)
U = pd.concat(U, ignore_index=True)

V = []
for i in range(nchan):
    s4 = pd.DataFrame([[np.average(data[0,3,i,ps:pf])]])
    V.append(s4)
V = pd.concat(V, ignore_index=True)


# calculate polarisation fraction
l = np.sqrt(np.float64(np.sqrt(Q**2+U**2)))/I
c = np.fabs(np.float64(V)/I)


x = np.arange(nchan)
xticks = np.linspace(0,3328, num=14).astype(int)
xticks_x = np.linspace(0,nchan, len(xticks))

plt.figure(figsize=(15,10),dpi=300)
ax = plt.subplot(111)
ax.set_ylim(0,1)
plt.scatter(x,l[0], label='Linear', c='r')
plt.scatter(x,c[0], label='Circular', c='b')
plt.xticks(xticks_x, xticks)
plt.xlabel('Frequency channel')
plt.ylabel('Polarisation Fraction')
plt.title('%s'%(sys.argv[1]))
plt.legend()

plt.savefig('F_PF_%s.png'%(sys.argv[1].split(os.extsep, 1)[0]))