import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import os

# python modulation_index.py <file> <phase start> <phase finish> <bin width>
# <bin width> is a power of 2

a = psrchive.Archive_load(sys.argv[1])
c = a.clone()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape

# on-pulse start and finish phase bins
on1 = float(sys.argv[2])
on2 = float(sys.argv[3])
on_s = int(on1*nbin)
on_f = int(on2*nbin)

b_width = int(sys.argv[4])
M = np.empty(0)
for i in np.arange(0,nbin,b_width):
    mean = np.average(data[0,0,0,i:(i+b_width)])
    rms = np.std(data[0,0,0,i:(i+b_width)])

    Mindex = rms/mean
    M = np.concatenate((M, Mindex), axis=None)

### PLOT R-FACTOR
ticks = np.round(np.linspace(on1,on2,num=11),2)
ticks_x = np.linspace(0,nbin/b_width-1,num=11)

plt.figure(figsize=(15,10),dpi=300)
ax = plt.subplot(111)
ax.set_xlim(0,nbin/b_width-1)
plt.plot(M, c='black')
#plt.axhline(25, color='red', linestyle='dashed', linewidth=1)
plt.xticks(ticks_x, ticks)
plt.xlabel('Phase')
plt.ylabel('RMS/Mean')
plt.title('%s, bin size = %s'%(sys.argv[1].split(os.extsep, 1)[0], b_width))

plt.savefig('modulation_index_%s_n%s.pdf'%(sys.argv[1].split(os.extsep, 1)[0], b_width))
print('modulation_index_%s_n%s.pdf'%(sys.argv[1].split(os.extsep, 1)[0], b_width))
