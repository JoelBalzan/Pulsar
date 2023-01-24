import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import psrchive
from scipy import odr
from scipy.signal import savgol_filter


a = psrchive.Archive_load(sys.argv[1])
c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape

p1 = float(sys.argv[2])
p2 = float(sys.argv[3])
ps = int(p1*nbin)
pf = int(p2*nbin)

x = np.arange(pf-ps)
y = data[0,0,0,ps:pf]
window = int(sys.argv[4])
window = int(np.ceil(window) // 2*2 +1)

y = savgol_filter(y, window, int(sys.argv[5]))

order = int(sys.argv[5])

poly_model = odr.polynomial(order)
d = odr.Data(x, y)
odr_obj = odr.ODR(d, poly_model)
output = odr_obj.run()
poly = np.poly1d(output.beta[::-1])
poly_y = poly(x)

ticks = np.round(np.linspace(p1,p2,num=11),4)
ticks_x = np.linspace(0,pf-ps+1,num=11)

plt.figure(figsize=(15,10),dpi=300)
ax = plt.subplot(111)
ax.set_xlim(0,pf-ps)
plt.xticks(ticks_x, ticks)
plt.plot(x, y, label="Flux Density", c='k')
plt.plot(x, poly_y, label="polynomial ODR")
plt.legend()
plt.title('%s, Polynomial Order = %s'%(sys.argv[1].split('.')[0], order))

plt.savefig('fit_poly_%s_n%s.pdf'%(sys.argv[1].split('.')[0], order))
print('fit_poly_%s_n%s.pdf'%(sys.argv[1].split('.')[0], order))