import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import psrchive
import glob


# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[1]

# Phase zoom factor
z = 0.0002


P = []
files = glob.glob("*.rescaled")
for ar in glob.glob("*.rescaled"):
    if p == "I":
        a = psrchive.Archive_load(ar)
        a = a.clone()
        a.remove_baseline()
        a.tscrunch()
        a.pscrunch()
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
    else:
        a = psrchive.Archive_load(ar)
        a = a.clone()
        a.remove_baseline()
        a.tscrunch()
        #a.bscrunch(8)
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

#### PLOTTING ####
x = int(len(files)*7)
fig = plt.figure(figsize=(x,10),dpi=300)
g = gridspec.GridSpec(ncols=len(files), nrows=1, wspace=0)

yticks = np.linspace(f1,f2, num=14).astype(int)
yticks_y = np.linspace(0,ff-fs-1, len(yticks))

for i in range(len(files)):
    ax = fig.add_subplot(g[0,i])
    ax.imshow(P[i],aspect='auto',cmap='viridis',origin='lower', 
              vmin=np.min(P[i]), 
	          vmax=0.5*np.max(P[i]), 
              interpolation='kaiser')
    ax.set_xticks([])
    ax.set_yticks(yticks_y)
    ax.set_title(files[i].split('.')[0],fontsize=10)

    if i == 0:
        ax.set_ylabel('Frequency (MHz)')
        ax.set_yticks(yticks_y)
        ax.set_yticklabels(yticks)
    else:
        ax.tick_params(bottom=False, labelbottom=False, left=True, labelleft=False, right=True, labelright=False)


### SAVE FIGURE
plt.savefig('peak_compare_PX500_38329.pdf', bbox_inches='tight')
print('peak_compare_PX500_38329.pdf')
