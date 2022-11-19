import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import psrfits_archive
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator

def vonMises (x, h, c, mu):
    # x: phase; h: height; c: concentration; mu: position
    return h*np.exp(c*(np.cos(2*np.pi*(x-mu))-1.0))


###########################################################

filename = sys.argv[1]  # file name of the fits file

#fig = plt.figure(figsize=(10,8))
#ax = fig.add_subplot(211)
#minorLocator = MultipleLocator(0.02)

fig = plt.figure(figsize=(12, 8)) 
ax0 = plt.subplot(111)
#gs = gridspec.GridSpec(3, 1, hspace=0, wspace=0.1, left=0.17) 
#ax0 = plt.subplot(gs[1:3])
#ax1 = plt.subplot(gs[0])

###########################################################

ar = psrfits_archive.archive(filename)
ar.get_profile(wts=False)
#print(ar.dat_scl, ar.dat_offs, ar.dat_wts, ar.profile[0,0,0,:])
#print(ar.dat_scl, ar.dat_offs, ar.dat_wts)
ar.remove_baseline(sub=0, delta=0.5)
print(ar.profile.shape)
ar.centering(sub=0)
#ar.fscrunch(16)
#ar.bscrunch(2)

ar.cal_pa(0, 0, delta=0.5)

#profile = np.rot90(np.flip(ar.profile))
#profile = np.rot90(np.flip(ar.profile))
profile = ar.profile
print(ar.profile.shape)

labeltxt = "{0:0.2f}".format(ar.mjd)

#print(ar.bw, ar.freq, ar.nchn)
#print(np.amin(profile), np.amax(profile))

################
minorLocator = MultipleLocator(0.02)
ax0.xaxis.set_minor_locator(minorLocator)

minorLocator = AutoMinorLocator(2)
ax0.yaxis.set_minor_locator(minorLocator)

ax0.set_ylim(-0.1, 1.5)
ax0.set_yticks(np.arange(0.0,1.5,0.2))
ax0.set_yticklabels(np.round(np.arange(0.0,1.5,0.2),1), fontsize=16)
ax0.set_xlim(0.0, 1)
ax0.set_xticks(np.arange(0,1.1,0.1))
ax0.set_xticklabels(np.round(np.linspace(0,1,11),1), fontsize=16)
#stokeI = np.squeeze(np.mean(profile, axis=0))
stokesI = profile[0,0,0,:]
phase = np.linspace(0,1,ar.nbin)
y = 0*phase

ax0.plot(phase, stokesI, ls='-', color='k', lw=3, label='Total intensity')
ax0.plot(phase, y, ls='--', color='k')

linear = np.mean(ar.linear, axis=0)
ax0.plot(phase, linear, ls='--', color='r', lw=3, label='Linear')
stokesV = profile[0,-1,0,:]
ax0.plot(phase, stokesV, ls='-.', color='b', lw=3, label='Circular')

ax0.set_xlabel('Pulse phase', fontsize=16)
ax0.set_ylabel('Flux density (mJy)', fontsize=16)
ax0.legend(loc='upper left', numpoints=1, fontsize=14)

#print (np.amax(stokesI)/1000., np.amin(stokesI)/1000.)
################
'''
ax1.tick_params(axis="x", which='both', direction="in")
minorLocator = MultipleLocator(0.02)
ax1.xaxis.set_minor_locator(minorLocator)
minorLocator = AutoMinorLocator(5)
ax1.yaxis.set_minor_locator(minorLocator)
#ax1.yaxis.set_minor_locator(minorLocator)
ax1.set_ylabel('Position angle (deg)', fontsize=16)
#ax1.set_ylim(-60, 10)
#ax1.set_yticks(np.arange(-45,10.,20))
#ax1.set_yticklabels(np.arange(-45,10.,20), fontsize=16)
ax1.set_xlim(0.0, 1)
ax1.set_xticks(np.arange(0,1.1,0.1))
ax1.set_xticklabels([])
ax1.errorbar(ar.phase_pa/ar.nbin, ar.pa, yerr=ar.err, fmt='.', color='k')
'''

###############
plt.savefig('profile.png', bbox_inches='tight')
#plt.show()
