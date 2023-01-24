import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import psrchive
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_widths


def func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y
    

a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.fscrunch()
a.pscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

p1 = float(sys.argv[2])
p2 = float(sys.argv[3])
ps = int(p1*nbin)
pf = int(p2*nbin)


# flux in Jy
flux = data[0,0,0,ps:pf]/1000
peaks, _ = find_peaks(flux)
fluxes = np.sort(flux[peaks])[::-1][0:2]

# index of highest 10 fluxes
idxs = []
for i in fluxes:
    idx = np.where(flux==i)[0][0]
    idxs.append(idx)
fluxes_i = np.array(idxs)

widths = peak_widths(flux, fluxes_i, rel_height=0.95)[0]

guess = []
for i in range(len(fluxes)):
    guess.append(fluxes_i[i])
    guess.append(fluxes[i])
    guess.append(widths[i])


popt, pcov = curve_fit(func, flux, np.arange(0,pf-ps,1), p0=guess)
fit = func(flux, *popt)


plt.figure(figsize=(15,10),dpi=300)
plt.plot(flux)
plt.plot(flux, fit, 'r-')
plt.plot(fluxes_i, fluxes, "x")
plt.savefig("fit_gaussian.pdf")


#ticks = np.round(np.linspace(p1,p2,num=11),2)
#ticks_x = np.linspace(0,pf-ps+1,num=11)
#
#plt.figure(figsize=(15,10),dpi=300)
#ax = plt.subplot(111)
#ax.set_xlim(0,pf-ps)
#
#
#plt.plot(data[0,0,0,ps:pf], label='Flux Density', c='black')
#plt.xticks(ticks_x, ticks)
#plt.xlabel('Phase')
#plt.ylabel('Flux Density')
#plt.title('%s'%sys.argv[1].split(os.extsep, 1)[0])
#plt.legend()
#
#plt.savefig("flux_pol_%s.pdf"%sys.argv[1].split(os.extsep, 1)[0])
#print("flux_pol_%s.pdf"%sys.argv[1].split(os.extsep, 1)[0])