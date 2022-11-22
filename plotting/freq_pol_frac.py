import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import pandas as pd

a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

p1 = float(sys.argv[2])
p2 = float(sys.argv[3])
ps = int(p1*nbin)
pf = int(p2*nbin)

#print(np.array([1,np.average(data[0,0,1,ps:pf])]))
I = []
for i in range(nchan):
    s1 = pd.DataFrame([[np.average(data[0,0,i,ps:pf])]])
    I.append(s1)
I = pd.concat(I, ignore_index=True)
I[I==0] = np.nan

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

L = []
for i in range(nchan):
    s23 = pd.DataFrame([[np.average(np.sqrt(data[0,1,i,ps:pf]**2+data[0,2,i,ps:pf]**2))]])
    L.append(s23)
L = pd.concat(L, ignore_index=True)

V = []
for i in range(nchan):
    s4 = pd.DataFrame([[np.average(data[0,3,i,ps:pf])]])
    V.append(s4)
V = pd.concat(V, ignore_index=True)


l = np.sqrt(np.float64(Q*Q+U*U))/I
l2 = np.sqrt(np.float64(L))/I
c = np.fabs(np.float64(V)/I)

x = np.arange(nchan)
plt.figure(figsize=(15,10),dpi=300)
plt.plot(x,l[0], label='l')
plt.plot(x,l2[0], label='l2')
plt.xlabel('Frequency channel')
plt.ylabel('Polarisation Fraction')
plt.title('Linear Polarisation')
plt.legend()

plt.savefig('F_PF.png')