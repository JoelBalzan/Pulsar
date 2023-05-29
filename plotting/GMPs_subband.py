import gc
import psrchive
import os
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import matplotlib.cm as cm

#project code
PCODE = sys.argv[1]
files = sorted(glob.glob("*.rescaled"))
if os.path.isfile('dynamic_spectra_'+PCODE+'.npy'):
    dynamic_spectra = np.load('dynamic_spectra_'+PCODE+'.npy')
else:
    dynamic_spectra = []
    counter = 0
    for ar in files:
        a = psrchive.Archive_load(ar)
        a.remove_baseline()
        a.tscrunch()
        a.fscrunch(8)
        a.pscrunch()
        a.centre()
        data = a.get_data()
        dynamic_spectra.append(data[0,0,:,:]/1000)

        # file progress counter
        counter += 1
        print("%s/%s"%(counter,len(files))," files completed", end='\r')
    ## SAVE FLUXES TO FILE
    np.save('dynamic_spectra_'+PCODE+'.npy', dynamic_spectra)
nfile, nchan, nbin = np.shape(dynamic_spectra)

### PLOTTING ###
# load first file to get period
a = psrchive.Archive_load(files[0])
# s per bin
period = a.integration_length()
spb = period/nbin

# mean dynamic spectrum
mean_DS = np.mean(dynamic_spectra, axis=0)

n_subbands = 26
n_chan_profile = int(nchan/n_subbands)

# brightest bin in each profile
brightest = []
mean_profiles = []
brightest_mean = []
GMP_subband = []
for i in range(n_subbands):
    brightest.append(np.amax(np.amax(dynamic_spectra[:,i*n_chan_profile:(i+1)*n_chan_profile,:], axis=0), axis=0))
    mean_profiles.append(np.mean(mean_DS[i*n_chan_profile:(i+1)*n_chan_profile,:], axis=0))
    brightest_mean.append(brightest[i]/np.max(mean_profiles[i]))
    if np.max(brightest_mean[i]) > 50:
        GMP_subband.append(i)
        print('GMP detected in subband %s'%(i+1))

del dynamic_spectra
del mean_DS
del brightest
gc.collect()

# phase start and finish
ps = 0#int(np.argmax(mean_profile)-0.12*nbin)
pf = 2**15#int(np.argmax(mean_profile)+0.12*nbin)
# trim arrays
#brightest_mean = brightest_mean[ps:pf]
#mean_profile = mean_profiles[ps:pf]

# xticks
xtickslabels = np.round(np.linspace(ps*spb, pf*spb, 11),3)
xticks = np.linspace(0, pf-ps-1, len(xtickslabels))

A4x, A4y = 8.27, 11.7
fig = plt.figure(figsize=(A4x,A4y),dpi=100)
g = gridspec.GridSpec(ncols=1, nrows=n_subbands, hspace=0)

for i in range(n_subbands):
    ax = fig.add_subplot(g[-1-i])
    scale = np.max(brightest_mean[i])/np.max(mean_profiles[i])
    ax.plot(mean_profiles[i]*scale, color='black', linewidth=0.5)
    ax.plot(brightest_mean[i], color='blue', linewidth=0.5, linestyle='--')
    ax.text(0.95, 0.95, i+1, transform=ax.transAxes, fontsize=10, fontweight='bold', verticalalignment='top', horizontalalignment='right', color='k')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtickslabels, fontsize=10)
    ax.tick_params(axis='both', which='major', direction='in', labelbottom=False, labeltop=False)
    ax.margins(x=0)
    if i == n_subbands-1:
        ax.tick_params(axis='both', which='major', direction='in', labelbottom=True, labeltop=False)
        ax.set_xlabel('Phase (s)', fontsize=8)
plt.savefig('GMP_subband_'+PCODE+'.pdf', dpi=100, bbox_inches='tight')
print('GMP_subband_'+PCODE+'.pdf')