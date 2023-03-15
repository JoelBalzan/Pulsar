import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import psrchive
from scipy.signal import find_peaks, peak_widths
import matplotlib.cm as cm


a = psrchive.Archive_load(sys.argv[1])
# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[2]

c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape


### FREQ ZOOM
f_scr = (4032-704)/a.get_nchan()
f1 = float(sys.argv[3])
f2 = float(sys.argv[4])
fs = int((f1-704)/f_scr)
ff = int((f2-704)/f_scr)

# neutron star period (ms)
period = a.integration_length()

# flux in Jy
flux = data[0,0,0,:]/1000
# minimum height of peaks (Jy)
h=3
peaks, _ = find_peaks(flux, height=h, distance=6)
# make sure there are an even number of peaks
if len(peaks)%2==1:
    # remove the lowest peak
    peaks = np.delete(peaks, np.argmin(flux[peaks]))

# width of peaks for setting imshow widths
w = np.round(peak_widths(flux, peaks, rel_height=0.8)[0]).astype(int)
for i in range(len(peaks)):
    if w[i] <= 4:
        w[i] = 5


# extract dynamic spectra and pulse profiles of sub-pulses
P = []
F = []
if p == "I":
    a.remove_baseline()
    a.tscrunch()
    a.pscrunch()
    data1 = a.get_data()
    nsub, npol, nchan, nbin = data1.shape

    for i in range(len(peaks)):    
        # total intensity
        I = data1[0,0,fs:ff,(peaks[i]-w[i]):(peaks[i]+w[i])]
        P.append(I)
        f = np.mean(I, axis=0)
        F.append(f)
else:
    a.remove_baseline()
    a.tscrunch()
    data1 = a.get_data()
    nsub, npol, nchan, nbin = data1.shape

    for i in range(len(peaks)):    
        if p == "SI":
            SI = data1[0,0,fs:ff,(peaks[i]-w[i]):(peaks[i]+w[i])]
            P.append(SI)
            f = np.mean(SI, axis=0)
            F.append(f)

        if p == "SQ":
            SQ = data1[0,1,fs:ff,(peaks[i]-w[i]):(peaks[i]+w[i])]
            P.append(SQ)
            f = np.mean(SQ, axis=0)
            F.append(f)

        if p == "SU":
            SU = data1[0,2,fs:ff,(peaks[i]-w[i]):(peaks[i]+w[i])]
            P.append(SU)
            f = np.mean(SU, axis=0)
            F.append(f)

        if p == "L":
            L = np.sqrt(data1[0,1,fs:ff,(peaks[i]-w[i]):(peaks[i]+w[i])]**2
                        +data1[0,2,fs:ff,(peaks[i]-w[i]):(peaks[i]+w[i])]**2)
            P.append(L)
            f = np.mean(L, axis=0)
            F.append(f)
            
        if p == "SV":
            SV = data1[0,3,fs:ff,(peaks[i]-w[i]):(peaks[i]+w[i])]
            P.append(SV)
            f = np.mean(SV, axis=0)
            F.append(f)


#### PLOTTING ####
x = int(len(peaks)*2)
fig = plt.figure(figsize=(x,30),dpi=300)
half_len = len(peaks)/2
g = gridspec.GridSpec(ncols=half_len, nrows=5, wspace=0, hspace=0, 
                      height_ratios=[2,1,7,1,7])

yticks = np.linspace(f1,f2, num=14).astype(int)
yticks_y = np.linspace(0,ff-fs-1, len(yticks))


# whole pulse profile
colours = cm.tab20(np.linspace(0, 1, len(peaks)))
ax0 = fig.add_subplot(g[0,:])
ax0.plot(flux,c='k')
ax0.set_xlim(left=peaks[0]-50, right=peaks[-1]+50)
ax0.set_ylabel('Flux (Jy)', fontsize=20)
ax0.tick_params(bottom=False, labelbottom=False, left=True, labelleft=True, labelsize=20)

# use the variance of the dynamic spectra to set the colour limits
var = []
for i in P:
    var.append(np.var(i))
var = var/np.min(var)

vmin = []
vmax = []
peak_plot = np.argmax(flux[peaks])
#min_plot = np.argmin(flux[peaks])
for i in var:
    if i < 1.5:
        vmin.append(0.15*np.min(P[peak_plot]))
        vmax.append(0.175*np.max(P[peak_plot]))
    else:
        vmin.append(0.3*np.min(P[peak_plot]))
        vmax.append(0.35*np.max(P[peak_plot]))

# plot dynamic spectra and sub-pulse profiles
for i in range(len(peaks)):
    # top flux density plot peaks
    ax0.plot(peaks[i], flux[peaks[i]], "*", ms=15., mew=0.8, c=colours[i])
    if i < half_len:
        # individual flux density plots row 1
        ax = fig.add_subplot(g[1,i])
        # sub-pulse profile 
        ax.plot(F[i],c='k')
        # sub-pulse profile peak markers
        pk = np.round(np.shape(F[i])[0]/2).astype(int)
        ax.plot(pk, F[i][pk], "*", ms=15., mew=0.8, c=colours[i])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(left=0, right=len(F[i])-1)
        ax.set_ylim(top=1.2*np.max(F[i]))

        # dynamic spectra row 1
        ax = fig.add_subplot(g[2,i])
        ax.imshow(P[i],aspect='auto',cmap='Spectral',origin='lower', interpolation='kaiser'
                  ,vmin=vmin[i], vmax=vmax[i]) 
        ax.set_xticks([])
        ax.set_yticks(yticks_y)

        if i == 0:
            ax.set_ylabel('Frequency (MHz)', fontsize=20)
            ax.set_yticks(yticks_y)
            ax.set_yticklabels(yticks, fontsize=20)
            ax.tick_params(direction='inout')
        else:
            ax.tick_params(bottom=False, labelbottom=False, left=True, labelleft=False, 
                           right=True, labelright=False, direction='inout')

    if i >= half_len:
        # individual flux density plots row 2
        ax = fig.add_subplot(g[3,i-half_len])
        # sub-pulse profile 
        ax.plot(F[i],c='k')
        # sub-pulse profile peak markers
        pk = np.round(np.shape(F[i])[0]/2).astype(int)
        ax.plot(pk, F[i][pk], "*", ms=15., mew=0.8, c=colours[i])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(left=0, right=len(F[i])-1)
        ax.set_ylim(top=1.2*np.max(F[i]))

        # dynamic spectra row 2
        ax = fig.add_subplot(g[4,i-half_len])
        ax.imshow(P[i],aspect='auto',cmap='Spectral',origin='lower', interpolation='kaiser'
                  ,vmin=vmin[i], vmax=vmax[i])

        ax.set_xticks([])
        ax.set_yticks(yticks_y)
        if i-half_len == 0:
            ax.set_ylabel('Frequency (MHz)', fontsize=20)
            ax.set_yticks(yticks_y)
            ax.set_yticklabels(yticks, fontsize=20)
            ax.tick_params(direction='inout')
        else:
            ax.tick_params(bottom=False, labelbottom=False, left=True, labelleft=False, 
                           right=True, labelright=False, direction='inout')





### SAVE FIGURE
plt.savefig('sub-pulse_compare_%s_%s.pdf'%(p,sys.argv[1].split('.')[0]), 
            bbox_inches='tight', dpi=300)
print('sub-pulse_compare_%s_%s.pdf'%(p,sys.argv[1].split('.')[0]))