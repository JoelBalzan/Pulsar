import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import psrchive
import os
from scipy.signal import find_peaks, peak_widths

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

a = psrchive.Archive_load(sys.argv[1])
# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[2]


### DETERMINE PEAK FLUX AND INDEX FOR PLOT CENTRING
c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape

# neutron star period (ms)
period = a.integration_length()


# peak and index
flux = data[0,0,0,:]
peaks, _ = find_peaks(flux)
peak_flux = np.sort(flux[peaks])[-1]
peak_idx = np.where(flux==peak_flux)[0][0]

# on-pulse phase start and finish
p1 = np.round(peak_idx/nbin - 0.001, 4)
p2 = np.round(peak_idx/nbin + 0.001, 4)


### PEAK INDICES  
# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))

flux = data[0,0,0,ps:pf]/1000
peaks, _ = find_peaks(flux)

# highest fluxes
fluxes = np.sort(flux[peaks])[::-1][0:6]

# index of highest fluxes
fluxes_i = []
for i in fluxes:
    idx = np.where(flux==i)[0][0]
    fluxes_i.append(idx)
fluxes_i = np.array(fluxes_i)

# peak widths
widths = peak_widths(flux, fluxes_i, rel_height=0.4)

# peak minimas
mins, _ = find_peaks(-flux)
# associate peaks with minimas
peak_mins = []
for i in fluxes_i:
    for j in range(len(mins)):
        if mins[j] < i < mins[j+1]:
            mins_i = np.array([[mins[j], mins[j+1]]])[0]
            peak_mins.append(mins_i)
peak_mins = np.array(peak_mins)


if sys.argv[2] == "I":
    c1 = a.clone()
    c1.remove_baseline()
    c1.tscrunch()
    c1.pscrunch()
    #c1.bscrunch(8)
    data2 = c1.get_data()
    nsub, npol, nchan, nbin = data2.shape

    # on-pulse phase bin start and finish
    ps = int(np.round(p1*nbin))
    pf = int(np.round(p2*nbin))

    # intensity
    I = data2[0,0,:,ps:pf]

    ### FREQ PEAKS
    freq_i = []
    for i in fluxes_i:
        idx = np.argmax(data2[0,0,:,i])
        freq_i.append(idx)
    freq_i = np.array(freq_i)

else:
    c1 = a.clone()
    c1.remove_baseline()
    c1.tscrunch()
    c1.bscrunch(8)
    data2 = c1.get_data()
    nsub, npol, nchan, nbin = data2.shape

    # on-pulse phase bin start and finish
    ps = int(np.round(p1*nbin))
    pf = int(np.round(p2*nbin))

    # polarisations
    SI = data2[0,0,:,ps:pf]
    SQ = data2[0,1,:,ps:pf]
    SU = data2[0,2,:,ps:pf]
    L = np.sqrt(data2[0,1,:,ps:pf]**2+data2[0,2,:,ps:pf]**2)
    SV = data2[0,3,:,ps:pf]

    ### FREQ PEAKS
    freq_i = []
    for i in fluxes_i:
        idx = np.argmax(eval(p)[:,i])
        freq_i.append(idx)
    freq_i = np.array(freq_i)


#### PLOTTING ####
fig = plt.figure(figsize=(10, 12), dpi=300) 
g = gridspec.GridSpec(5, 5, hspace=0, wspace=0, left=0.17) 
ax0 = plt.subplot(g[0,:4])
ax1 = plt.subplot(g[1:5, :4])
ax2 = plt.subplot(g[1:5, 4])
ax0.tick_params(axis='both', which='major', labelsize=12)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)

# seconds per bin
bs = 1000*period/nbin
nbin_zoom = np.shape(eval(p))[1]


### PLOT SPECTRUM
spectrum = np.mean(eval(p)[:,:], axis=1)
x = 0*np.arange(nchan)

ax2.set_ylim(0.0, nchan-1)
ax2.set_xlim(-3000,np.max(spectrum)+1000)
ax2.set_yticks(np.arange(0, nchan, 16))
ax2.set_yticklabels([])

n_ytick = 8
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
minorLocator = MultipleLocator(nchan/n_ytick/7.)
ax2.yaxis.set_minor_locator(minorLocator)
ax2.set_yticks(np.linspace(0.0, nchan, n_ytick))
freq0 = a.get_centre_frequency() - a.get_bandwidth()/2.
cfreq = a.get_bandwidth()/nchan
ax2.set_yticklabels([])
ax2.tick_params(axis="y", which='both', direction="in", pad=-22)

ax2.plot(spectrum, np.arange(nchan), ls='-', color='k', lw=2)
ax2.plot(x, np.arange(nchan), ls='--', color='k')


### PLOT POLARISATION
xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.linspace(704,4032, num=14).astype(int)
yticks_y = np.linspace(0,nchan, len(yticks))


### MASK ZAPPED CHANNELS
#masked_data = np.ma.masked_values(eval(p), 0.)
#cmap = matplotlib.cm.get_cmap("Spectral").copy()
#cmap.set_bad(color='white')
#mask zapped channels colour
#if len(sys.argv)==5:
#    pol = np.ma.masked_values(eval(p), 0.)
#    cmap = matplotlib.cm.get_cmap("Spectral").copy()
#    cmap.set_bad(color='white')
#else:
#    pol = eval(p)

# order peak indices
fluxes_i, freq_i = zip(*sorted(zip(fluxes_i, fluxes, freq_i)))

ax1.imshow(eval(p), cmap="Spectral", vmin=np.min(eval(p)), vmax=0.3*np.max(eval(p)), aspect='auto', origin='lower')
#ax1.vlines(widths[2][0], 0, nchan, colors='k', linestyles='dashed')
#ax1.vlines(widths[3][0], 0, nchan, colors='k', linestyles='dashed')
# plot peaks on dynamic spectra
#ax1.plot(fluxes_i, freq_i, marker='x', c='k')

ax1.set_xlim(0.0, pf-ps-1)
ax1.set_xticks(xticks_x)
ax1.set_xticklabels(xticks, fontsize=12)
ax1.set_xlabel('Time (ms)', fontsize=14)

ax1.set_ylabel('Frequency (MHz)', fontsize=12)
ax1.set_ylim(0, nchan-1)
ax1.set_yticks(yticks_y)
ax1.set_yticklabels(yticks, fontsize=12)


### PLOT FLUX DENSITY
nsub, npol, nchan, nbin = data.shape

# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))

xticks = np.round(np.linspace((-nbin_zoom/2)*bs,(nbin_zoom/2)*bs,num=11),2)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
yticks = np.round(np.linspace(0,peak_flux/1000 - 1, num=4)).astype(int)
yticks_y = np.linspace(0,(peak_flux/1000)-1,num=len(yticks))

ax0.plot(flux, c='k', lw=2)
# plot peaks and mins
#ax0.plot(fluxes_i, flux[fluxes_i], 'x', c='b')
#ax0.plot(mins, flux[mins], 'x', c='r')
ax0.set_xlim(0.0, pf-ps-1)
ax0.set_ylim(-5,peak_flux/1000+5)
ax0.tick_params(axis="x", which='both', direction="in", pad=-22)
ax0.set_xticks(xticks_x)
ax0.set_xticklabels([])
#ax0.set_yticks([])
phase = np.linspace(0,1,nbin)
y = 0*np.arange(nbin)
ax0.plot(np.arange(nbin), y, ls='--', color='k')
ax0.set_ylabel('Flux density (Jy)', fontsize=12)

ax0.set_title('%s Polarisation %s'%(p,sys.argv[1].split('.')[0]), fontsize=12)




### SAVE FIGURE
plt.savefig("pulse_drift_%s.pdf"%(sys.argv[1].split('.')[0]))
print("pulse_drift_%s.pdf"%(sys.argv[1].split('.')[0]))


