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
        a.fscrunch(4)
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
print(nfile, nchan, nbin)

### PLOTTING ###
# load first file to get period
a = psrchive.Archive_load(files[0])
# s per bin
period = a.integration_length()
spb = period/nbin

# mean dynamic spectrum
mean_DS = np.mean(dynamic_spectra, axis=0)
print(np.shape(mean_DS))

n_subbands = 26
n_chan_profile = nchan/n_subbands

# brightest bin in each profile
brightest = []
for i in range(nfile):
    for j in range(n_subbands):
        brightest.append(np.amax(dynamic_spectra[i,j*n_chan_profile:(j+1)*n_chan_profile,:], axis=0))

mean_profiles = []
for i in range(n_subbands):
    mean_profiles.append(np.mean(mean_DS[i*n_chan_profile:(i+1)*n_chan_profile,:], axis=0))

brightest_mean = brightest/mean_profiles

# phase start and finish
ps = 0#int(np.argmax(mean_profile)-0.12*nbin)
pf = 2**15#int(np.argmax(mean_profile)+0.12*nbin)
# trim arrays
brightest_mean = brightest_mean[ps:pf]
mean_profile = mean_profiles[ps:pf]

# xticks
xtickslabels = np.round(np.linspace(ps*spb, pf*spb, 11),3)
xticks = np.linspace(0, pf-ps-1, len(xtickslabels))

A4x, A4y = 8.27, 11.7
fig = plt.figure(figsize=(A4x,A4y),dpi=300)
g = gridspec.GridSpec(ncols=1, nrows=n_subbands, hspace=0)

for i in range(n_subbands):
    ax = fig.add_subplot(g[i])
    scale = np.max(brightest_mean[i])/np.max(mean_profiles[i])
    ax.plot(mean_profiles[i]*scale, color='black', linewidth=0.5)
    ax.plot(brightest_mean[i], color='blue', linewidth=0.5, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=8, direction='in', labelbottom=False, labeltop=False)
    if i == n_subbands-1:
        ax.tick_params(axis='both', which='major', labelsize=8, direction='in', labelbottom=True, labeltop=False)
        ax.set_xlabel('Phase (s)', fontsize=8)