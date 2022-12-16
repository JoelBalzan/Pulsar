import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
import psrchive
import os
import glob
from scipy.signal import find_peaks

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)


# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[1]
# neutron star period (ms)
period = 1/float(sys.argv[2])


#fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(15,15), dpi=300)

P = []
for ar in glob.glob("*.rescaled"):
    if sys.argv[1] == "I":
        a = psrchive.Archive_load(ar)
       
        c = a.clone()
        c.remove_baseline()
        c.tscrunch()
        c.fscrunch()
        c.pscrunch()
        data = c.get_data()
        nsub, npol, nchan, nbin = data.shape

        # peak and index
        flux = data[0,0,0,:]
        peaks, _ = find_peaks(flux)
        peak_flux = np.sort(flux[peaks])[-1]
        peak_idx = np.where(flux==peak_flux)[0][0]

        # on-pulse phase start and finish
        p1 = np.round(peak_idx/nbin - 0.006, 4)
        p2 = np.round(peak_idx/nbin + 0.006, 4)

        ### POLARISATION
        c1 = a.clone()
        c1.remove_baseline()
        c1.tscrunch()
        c1.pscrunch()
        data1 = c1.get_data()
        nsub, npol, nchan, nbin = data1.shape

        # on-pulse phase bin start and finish
        ps = int(np.round(p1*nbin))
        pf = int(np.round(p2*nbin))

        # intensity
        I = data1[0,0,:,ps:pf]

        P.append(I)

    else:
        a = psrchive.Archive_load(ar)
       
        c = a.clone()
        c.remove_baseline()
        c.tscrunch()
        c.fscrunch()
        c.pscrunch()
        data = c.get_data()
        nsub, npol, nchan, nbin = data.shape

        # peak and index
        flux = data[0,0,0,:]
        peaks, _ = find_peaks(flux)
        peak_flux = np.sort(flux[peaks])[-1]
        peak_idx = np.where(flux==peak_flux)[0][0]

        # on-pulse phase start and finish
        p1 = np.round(peak_idx/nbin - 0.006, 4)
        p2 = np.round(peak_idx/nbin + 0.006, 4)

        ### POLARISATION
        c1 = a.clone()
        c1.remove_baseline()
        c1.tscrunch()
        data1 = c1.get_data()
        nsub, npol, nchan, nbin = data1.shape

        # on-pulse phase bin start and finish
        ps = int(np.round(p1*nbin))
        pf = int(np.round(p2*nbin))

        # polarisations
        SI = data1[0,0,:,ps:pf]
        SQ = data1[0,1,:,ps:pf]
        SU = data1[0,2,:,ps:pf]
        L = np.sqrt(data1[0,1,:,ps:pf]**2+data1[0,2,:,ps:pf]**2)
        SV = data1[0,3,:,ps:pf]

        P.append(eval(p))

#### PLOTTING ####
fig = plt.figure(figsize=(15,15),dpi=300)
g = gridspec.GridSpec(ncols=3, nrows=3, hspace=0.08, wspace=0.08)

# milliseconds per bin
bs = 1000*period/nbin
nbin_zoom = np.shape(P[0])[1]

### PLOT POLARISATION
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=5)).astype(int)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.linspace(704,4032, num=14).astype(int)
yticks_y = np.linspace(0,nchan-1, len(yticks))


for i in range(9):
    ax = fig.add_subplot(g[i])
    ax.imshow(P[i], cmap="Spectral", vmin=np.min(P[i]), vmax=np.max(P[i]), aspect='auto', origin='lower')
    plt.yticks(yticks_y, yticks, fontsize=10)
    plt.xticks(xticks_x, xticks, fontsize=10)
    if i == 0:
        ax.set_ylabel('Frequency (MHz)', fontsize=10)
        ax.tick_params(bottom=True, labelbottom=False, left=True, labelleft=True, right=True, top=True)
    if i == 1:
        ax.tick_params(bottom=True, labelbottom=False, left=True, labelleft=False, right=True, top=True)
    if i == 2:
        ax.tick_params(bottom=True, labelbottom=False, left=True, labelleft=False, right=True, top=True)
    if i == 3:
        ax.set_ylabel('Frequency (MHz)', fontsize=10)
        ax.tick_params(bottom=True, labelbottom=False, left=True, labelleft=True, right=True, top=True)
    if i == 4:
        ax.tick_params(bottom=True, labelbottom=False, left=True, labelleft=False, right=True, top=True)
    if i == 5:
        ax.tick_params(bottom=True, labelbottom=False, left=True, labelleft=False, right=True, top=True)
    if i == 6:
        ax.set_xlabel('Time (ms)', fontsize=10)
        ax.set_ylabel('Frequency (MHz)', fontsize=10)
        ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True, right=True, top=True)
    if i == 7:
        ax.set_xlabel('Time (ms)', fontsize=10)
        ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=False, right=True, top=True)
    if i == 8:
        ax.set_xlabel('Time (ms)', fontsize=10)
        ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=False, right=True, top=True)


### SAVE FIGURE
plt.savefig('PWZ_%s_PX500_38329.pdf'%p)
print('PWZ_%s_PX500_38329.pdf'%p)
