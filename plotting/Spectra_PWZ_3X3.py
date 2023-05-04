import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import psrchive
import os
import glob
import string

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[1]

# Phase zoom factor
z = 0.0002


S = []
files = sorted(glob.glob("*.rescaled"))
# exit if <12 files
if len(files) < 12:
    print("12 files required. You only have %s files."%len(files))
    sys.exit()

for ar in files[0:len(files)]:
    if sys.argv[1] == "I":
        a = psrchive.Archive_load(ar)
        a.remove_baseline()
        a.tscrunch()
        a.pscrunch()
        #c1.bscrunch(8)
        data1 = a.get_data()
        nsub, npol, nchan, nbin = data1.shape

        # peak and index
        peak_idx = np.argmax(np.mean(data1[0,0,:,:], axis=0))

        # on-pulse phase start and finish
        p1 = np.round(peak_idx/nbin - z, 4)
        p2 = np.round(peak_idx/nbin + z, 4)

        # on-pulse phase bin start and finish
        ps = int(np.round(p1*nbin))
        pf = int(np.round(p2*nbin))

        ### FREQ ZOOM
        f_scr = (4032-704)/a.get_nchan()

        f1 = float(sys.argv[2])
        f2 = float(sys.argv[3])
        fs = int((f1-704)/f_scr)
        ff = int((f2-704)/f_scr)

        # intensity
        I = np.mean(data1[0,0,fs:ff,ps:pf]/1000, axis=1)

        S.append(I)

    else:
        a = psrchive.Archive_load(ar)
        a.remove_baseline()
        a.tscrunch()
        a.bscrunch(8)
        data1 = a.get_data()
        nsub, npol, nchan, nbin = data1.shape

        # peak and index
        peak_idx = np.argmax(data1.mean(axis=(1,2))[0])

        # on-pulse phase start and finish
        p1 = np.round(peak_idx/nbin - z, 4)
        p2 = np.round(peak_idx/nbin + z, 4)

        # on-pulse phase bin start and finish
        ps = int(np.round(p1*nbin))
        pf = int(np.round(p2*nbin))

        ### FREQ ZOOM
        f_scr = (4032-704)/a.get_nchan()

        f1 = float(sys.argv[2])
        f2 = float(sys.argv[3])
        fs = int((f1-704)/f_scr)
        ff = int((f2-704)/f_scr)

        # polarisations
        if p == "SI":
            SI = np.mean(data1[0,0,fs:ff,ps:pf]/1000, axis=1)
            S.append(SI)
        if p == "SQ":
            SQ = np.mean(data1[0,1,fs:ff,ps:pf]/1000, axis=1)
            S.append(SQ)
        if p == "SU":
            SU = data1[0,2,fs:ff,ps:pf]
            S.append(SU)
        if p == "L":
            L = np.mean(np.sqrt((data1[0,1,fs:ff,ps:pf]/1000)**2+(data1[0,2,fs:ff,ps:pf]/1000)**2), axis=1)
            S.append(L)
        if p == "SV":
            SV = np.mean(data1[0,3,fs:ff,ps:pf]/1000, axis=1)
            S.append(SV)


#### PLOTTING ####
bot = False

A4x = 8.27
A4y = 11.69
fig = plt.figure(figsize=(A4x,A4y),dpi=600)
g = gridspec.GridSpec(ncols=3, nrows=4, hspace=0.11, wspace=0.11)


### PLOT DYNAMIC SPECTRA ###
xticks = np.linspace(f1,f2, num=7).astype(int)
xticks_x = np.linspace(0,ff-fs-1, len(xticks))

# alphabet for plot labelling
alphabet = list(string.ascii_lowercase)

fontsize = 10
for i in range(len(files)):
    ax = fig.add_subplot(g[i])
    ax.plot(S[i], c='k', lw=1)
    ax.set_xticks(xticks_x[i])
    ax.set_xticklabels(xticks[i], fontsize=fontsize)
    if i == 0 or i == 3 or i == 6 or i == 9:
        ax.set_ylabel('Flux Density (Jy)', fontsize=fontsize)
        ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=True, right=True, top=True)
        ax.text(0.05, 0.95, alphabet[i]+")", transform=ax.transAxes, fontsize=fontsize, fontweight='bold', va='top', color='w')
    if i == 1 or i == 2 or i == 4 or i == 5 or i == 7 or i == 8:
        ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=False, right=True, top=True)
        ax.text(0.05, 0.95, alphabet[i]+")", transform=ax.transAxes, fontsize=fontsize, fontweight='bold', va='top', color='w')
    if i == 9:
        ax.set_xlabel('Frequency (MHz)', fontsize=fontsize)
        ax.set_ylabel('Flux Density (Jy)', fontsize=fontsize)
        ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True, right=True, top=True)
    if i == 10 or i == 11:
        ax.set_xlabel('Frequency (MHz)', fontsize=fontsize)
        ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=False, right=True, top=True)
        ax.text(0.05, 0.95, alphabet[i]+")", transform=ax.transAxes, fontsize=fontsize, fontweight='bold', va='top', color='w')


### SAVE FIGURE
plt.savefig('PWZ_%s_PX500.png'%p, bbox_inches='tight')
print('PWZ_%s_PX500.png'%p)
