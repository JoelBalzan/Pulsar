import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import psrchive
from scipy.signal import find_peaks, peak_widths
import matplotlib.cm as cm


### DEFINE FUNCTIONS ###
# plot dynamic spectra and sub-pulse profiles with 1 row for <= 5 peaks
def plot_subpulses_1row():
    x = int(len(peaks)*4)
    fig = plt.figure(figsize=(x,15),dpi=600)

    length = int(len(peaks))
    g = gridspec.GridSpec(ncols=length, nrows=3, wspace=0, hspace=0, 
                      height_ratios=[2,1,7])
    
    colours = cm.tab20(np.linspace(0, 1, len(peaks)))
    ax0 = fig.add_subplot(g[0,:])
    ax0.plot(flux, c='k', label=f"{np.round(mspb*(peaks[-1] - peaks[0] + 100), 2)} ms")
    ax0.set_xlim(left=peaks[0]-50, right=peaks[-1]+50)
    ax0.set_ylabel('Flux (Jy)', fontsize=30)
    ax0.tick_params(bottom=False, labelbottom=False, left=True, labelleft=True, labelsize=30)
    ax0.legend(fontsize=15, loc='best')
    for i in range(len(peaks)):
        # top flux density plot peaks
        ax0.plot(peaks[i], flux[peaks[i]], "*", ms=15., mew=0.8, c=colours[i])

        # individual flux density plots row 1
        ax = fig.add_subplot(g[1,i])
        # sub-pulse profile 
        ax.plot(F[i],c='k')
        # sub-pulse profile peak markers
        # centre channel
        pk = np.round(np.shape(F[i])[0]/2).astype(int)
        ax.plot(pk, F[i][pk], "*", ms=15., mew=0.8, c=colours[i])
        # plot peak widths
        ps, _ = find_peaks(F[i], height=h, distance=100)
        width = peak_widths(F[i], ps, rel_height=0.5)
        plt.hlines(*width[1:], color='r', lw=2, label = f"{np.round(mspb*width[0][0], 2)} ms")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(left=0, right=len(F[i])-1)
        ax.set_ylim(top=1.2*np.max(F[i]))
        ax.legend(fontsize=15, handlelength = 1, borderpad = 0.15, loc='best')

        # dynamic spectra row 1
        ax = fig.add_subplot(g[2,i])
        ax.imshow(P[i],aspect='auto',cmap='viridis',origin='lower', interpolation='kaiser'
                  ,vmin=vmin[i], vmax=vmax[i]) 
        ax.set_xticks([])
        ax.set_yticks(yticks_y)

        # ticks and labels
        if i == 0:
            ax.set_ylabel('Frequency (MHz)', fontsize=30)
            ax.set_yticks(yticks_y)
            ax.set_yticklabels(yticks, fontsize=30)
            ax.tick_params(direction='inout')
        else:
            ax.tick_params(bottom=False, labelbottom=False, left=True, labelleft=False, 
                           right=True, labelright=False, direction='inout')

# plot dynamic spectra and sub-pulse profiles with 2 rows for >5 peaks
def plot_subpulses_2row():
    x = int(len(peaks)*2)
    fig = plt.figure(figsize=(x,30),dpi=600)

    length = int(len(peaks)/2)
    g = gridspec.GridSpec(ncols=length, nrows=5, wspace=0, hspace=0, 
                      height_ratios=[2,1,7,1,7])
    
    colours = cm.tab20(np.linspace(0, 1, len(peaks)))
    ax0 = fig.add_subplot(g[0,:])
    ax0.plot(flux, c='k', label=f"{np.round(mspb*(peaks[-1] - peaks[0] + 100), 2)} ms")
    ax0.set_xlim(left=peaks[0]-50, right=peaks[-1]+50)
    ax0.set_ylabel('Flux (Jy)', fontsize=30)
    ax0.tick_params(bottom=False, labelbottom=False, left=True, labelleft=True, labelsize=30)
    ax0.legend(fontsize=15, loc='best')
    for i in range(len(peaks)):
        # top flux density plot peaks
        ax0.plot(peaks[i], flux[peaks[i]], "*", ms=15., mew=0.8, c=colours[i])
        if i < length:
            # individual flux density plots row 1
            ax = fig.add_subplot(g[1,i])
            # sub-pulse profile 
            ax.plot(F[i],c='k')
            # sub-pulse profile peak markers
            # centre channel
            pk = np.array([np.round(np.shape(F[i])[0]/2).astype(int)])
            ax.plot(pk, F[i][pk[0]], "*", ms=15., mew=0.8, c=colours[i])
            # plot peak widths
            width = peak_widths(F[i], pk, rel_height=0.5)
            plt.hlines(*width[1:], color='r', lw=2, label = f"{np.round(mspb*width[0][0], 2)} ms")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(left=0, right=len(F[i])-1)
            ax.set_ylim(top=1.2*np.max(F[i]))
            ax.legend(fontsize=15, handlelength = 1, borderpad = 0.15, loc='best')
    
            # dynamic spectra row 1
            ax = fig.add_subplot(g[2,i])
            ax.imshow(P[i],aspect='auto',cmap='viridis',origin='lower', interpolation='kaiser'
                      ,vmin=vmin[i], vmax=vmax[i]) 
            ax.set_xticks([])
            ax.set_yticks(yticks_y)
    
            if i == 0:
                ax.set_ylabel('Frequency (MHz)', fontsize=30)
                ax.set_yticks(yticks_y)
                ax.set_yticklabels(yticks, fontsize=30)
                ax.tick_params(direction='inout')
            else:
                ax.tick_params(bottom=False, labelbottom=False, left=True, labelleft=False, 
                               right=True, labelright=False, direction='inout')
        if i >= length:
            # individual flux density plots row 2
            ax = fig.add_subplot(g[3,i-length])
            # sub-pulse profile 
            ax.plot(F[i],c='k')
            # sub-pulse profile peak markers
            # centre channel
            pk = np.array([np.round(np.shape(F[i])[0]/2).astype(int)])
            ax.plot(pk, F[i][pk[0]], "*", ms=15., mew=0.8, c=colours[i])
            # plot peak widths
            width = peak_widths(F[i], pk, rel_height=0.5)
            plt.hlines(*width[1:], color='r', lw=2, label = f"{np.round(mspb*width[0][0], 2)} ms")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(left=0, right=len(F[i])-1)
            ax.set_ylim(top=1.2*np.max(F[i]))
            ax.legend(fontsize=15, handlelength = 1, borderpad = 0.15, loc='best')
    
            # dynamic spectra row 2
            ax = fig.add_subplot(g[4,i-length])
            ax.imshow(P[i],aspect='auto',cmap='viridis',origin='lower', interpolation='kaiser'
                      ,vmin=vmin[i], vmax=vmax[i])
    
            ax.set_xticks([])
            ax.set_yticks(yticks_y)
            if i-length == 0:
                ax.set_ylabel('Frequency (MHz)', fontsize=30)
                ax.set_yticks(yticks_y)
                ax.set_yticklabels(yticks, fontsize=30)
                ax.tick_params(direction='inout')
            else:
                ax.tick_params(bottom=False, labelbottom=False, left=True, labelleft=False, 
                               right=True, labelright=False, direction='inout')




### LOAD DATA ###
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

# neutron star period (s)
period = a.integration_length()
# seconds per bin
spb = period/nbin
# milliseconds per bin
mspb = spb*1000


# flux in Jy
flux = data[0,0,0,:]/1000
# minimum height of peaks (Jy)
h=3
peaks, _ = find_peaks(flux, height=h, distance=8)
# make sure there are an even number of peaks
if (len(peaks) != 1) and (len(peaks)%2 == 1) and (len(peaks) > 5):
    # remove the lowest peak
    peaks = np.delete(peaks, np.argmin(flux[peaks]))

# width of peaks for setting imshow widths
w = np.round(peak_widths(flux, peaks, rel_height=0.5)[0]).astype(int)
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
# use the variance of the dynamic spectra to set the colour limits
var = []
for i in P:
    var.append(np.var(i))
var = var/np.max(var)

vmin = []
vmax = []
peak_plot = np.argmax(flux[peaks])
for i in range(len(peaks)):
    vmin.append(0.8*var[i]*np.min(P[peak_plot]))
    vmax.append(0.5*var[i]*np.max(P[peak_plot]))


yticks = np.linspace(f1,f2, num=14).astype(int)
yticks_y = np.linspace(0,ff-fs-1, len(yticks))

# plot dynamic spectra
if len(peaks) <= 5:
    plot_subpulses_1row()
else:
    plot_subpulses_2row()


### SAVE FIGURE
plt.savefig('%s_%s_%s.pdf'%(sys.argv[0].split('.')[0], p, sys.argv[1].split('.')[0]), 
            bbox_inches='tight', dpi=600)
print('%s_%s_%s.pdf'%(sys.argv[0].split('.')[0], p, sys.argv[1].split('.')[0]))
