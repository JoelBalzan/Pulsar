import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import psrchive
import os
import glob

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)


# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[1]

# Phase zoom factor
z = 0.0002


P = []
xticks = []
xticks_x = []
files = sorted(glob.glob("*.pazi.rescaled"))
for ar in files[0:9]:
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
        I = data1[0,0,fs:ff,ps:pf]

        P.append(I)

        # milliseconds per bin
        period = a.integration_length()
        bs = period/nbin
        nbin_zoom = np.shape(eval(p))[1]

        xt = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=5),2)#.astype(int)
        xticks.append(xt)
        xt_x = np.linspace(0,pf-ps-1,num=len(xt))
        xticks_x.append(xt_x)

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
            SI = data1[0,0,fs:ff,ps:pf]
            P.append(SI)
        if p == "SQ":
            SQ = data1[0,1,fs:ff,ps:pf]
            P.append(SQ)
        if p == "SU":
            SU = data1[0,2,fs:ff,ps:pf]
            P.append(SU)
        if p == "L":
            L = np.sqrt(data1[0,1,fs:ff,ps:pf]**2+data1[0,2,fs:ff,ps:pf]**2)
            P.append(L)
        if p == "SV":
            SV = data1[0,3,fs:ff,ps:pf]
            P.append(SV)

        # milliseconds per bin
        period = a.integration_length()
        bs = period/nbin
        nbin_zoom = np.shape(eval(p))[1]

        xt = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=5),2)#.astype(int)
        xticks.append(xt)
        xt_x = np.linspace(0,pf-ps-1,num=len(xt))
        xticks_x.append(xt_x)

#### PLOTTING ####
fig = plt.figure(figsize=(15,15),dpi=300)
g = gridspec.GridSpec(ncols=3, nrows=3, hspace=0.15, wspace=0.15)


### PLOT POLARISATION
if np.all(xticks == xticks[0]):
    bot = False
else:
    bot = True

yticks = np.linspace(f1,f2, num=14).astype(int)
yticks_y = np.linspace(0,ff-fs-1, len(yticks))

var = []
for i in P:
    var.append(np.var(i))
    #print(10*np.std(np.mean(i, axis=0))/np.max(np.mean(i, axis=0)))
var = var/np.max(var)

vmin = []
vmax = []
for i in range(9):
    vmin.append(0.8*var[i]*np.min(P[i]))
    vmax.append(0.8*var[i]*np.max(P[i]))


for i in range(9):
    ax = fig.add_subplot(g[i])
    ax.imshow(P[i], cmap="viridis", vmin=vmin[i], vmax=vmax[i], aspect='auto', origin='lower', interpolation='kaiser')
    ax.set_xticks(xticks_x[i])
    ax.set_xticklabels(xticks[i], fontsize=12)
    plt.yticks(yticks_y, yticks, fontsize=12)
    if i == 0:
        ax.set_ylabel('Frequency (MHz)', fontsize=12)
        ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=True, right=True)
    if i == 1:
        ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=False, right=True)
    if i == 2:
        ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=False, right=True)
    if i == 3:
        ax.set_ylabel('Frequency (MHz)', fontsize=12)
        ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=True, right=True)
    if i == 4:
        ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=False, right=True)
    if i == 5:
        ax.tick_params(bottom=True, labelbottom=bot, left=True, labelleft=False, right=True)
    if i == 6:
        ax.set_xlabel('Time (s)', fontsize=12)
        ax.set_ylabel('Frequency (MHz)', fontsize=12)
        ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True, right=True)
    if i == 7:
        ax.set_xlabel('Time (s)', fontsize=12)
        ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=False, right=True)
    if i == 8:
        ax.set_xlabel('Time (s)', fontsize=12)
        ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=False, right=True)



### SAVE FIGURE
plt.savefig('PWZ_%s_PX500_38907.pdf'%p, bbox_inches='tight')
print('PWZ_%s_PX500_38907.pdf'%p)
