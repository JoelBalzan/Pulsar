import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import psrchive
import os
from scipy.signal import find_peaks, peak_widths
from scipy.signal import savgol_filter


# python pulse_drift.py <file> <Polarisation> <freq_window> <poly_degree>
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)


def cal_fwtm (freq, spec):
	peak = np.amax(spec)
	peak_idx = np.argmax(spec)
	tenth = 0.1*peak
	signs = np.sign(np.add(spec, -tenth))
	#print (peak, tenth, signs)
	zero_crossings = (signs[0:-1] != signs[1:])
	#zero_crossings = (signs[0:-2] != signs[1:-1])
	zero_crossings_i = np.where(zero_crossings)[0]
	
	#print (peak_idx)
	#print (zero_crossings, np.where(zero_crossings), zero_crossings_i, peak_idx)
	if len(zero_crossings_i) > 0:
		signs = np.sign(np.add(zero_crossings_i, -peak_idx))
		num = len(signs)
		#print (signs)
		if np.all(np.equal(signs, np.ones(num))):
			idx1 = 0
			idx2 = zero_crossings_i[0]
			#print ('here1')
		elif np.all(np.equal(signs, -np.ones(num))):
			idx1 = zero_crossings_i[-1]
			idx2 = -1
			#print ('here2')
		else:
			#print ('here3')
			#print (signs)
			#print (signs[0:-1], signs[1:])
			crossings2 = (signs[0:-1] != signs[1:])
			#crossings2 = (signs[0:-2] != signs[1:-1])
			#print (crossings2)
			crossings2_i = np.where(crossings2)[0]
			#print (crossings2_i)
		
			idx1 = zero_crossings_i[crossings2_i[0]]
			idx2 = zero_crossings_i[crossings2_i[0]+1]
	else:
		idx1 = 0
		idx2 = -1

	#print (idx1, idx2)
	freq1 = freq[idx1]
	freq2 = freq[idx2]
	#print (freq1, freq2)
	bw = np.fabs(freq1-freq2)

	#return freq1, freq2, np.fabs(freq1-freq2)
	return int(idx1), int(idx2), bw	


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
#p1 = np.round(peak_idx/nbin - 0.0008, 4)
#p2 = np.round(peak_idx/nbin + 0.0008, 4)
p1 = np.round(peak_idx/nbin - 0.001, 4)
p2 = np.round(peak_idx/nbin + 0.001, 4)

### PEAK INDICES  
# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))

flux = data[0,0,0,ps:pf]/1000
h=3
peaks, _ = find_peaks(flux, height=h)

# peak minimas
mins, _ = find_peaks(-flux, wlen=2)
# associate peaks with minimas
if peaks[-1] > mins[-1]:
	peaks = peaks[:-1]
if peaks[0] < mins[0]:
	peaks = peaks[1:]

peak_mins = []
for i in peaks:
	for j in range(len(mins)):
		if mins[j] < i < mins[j+1]:
			mins_i = np.array([[mins[j], mins[j+1]]])[0]
			peak_mins.append(mins_i)
peak_mins = np.array(peak_mins)


#### FIT POLYNOMIALS TO SPECTRA ####
### SPECTRA OF PEAKS
c1 = a.clone()
c1.remove_baseline()
c1.tscrunch()
c1.pscrunch()
data1 = c1.get_data()
nsub, npol, nchan, nbin = data1.shape
spectra = data1[0,0,:,ps:pf]/1000

### DEFINE SPECTRA VARIABLES
dict={}
for i in range(len(peaks)):
	key = str("S"+str(i))
	spec = np.mean(spectra[:,peak_mins[i][0]:peak_mins[i][1]], axis=1)
	dict[key] = spec.tolist()
for key,value in dict.items():
	exec(f'{key}={value}')


#### PLOTTING ####
freq_window = float(sys.argv[3])
freq_poly = int(sys.argv[4])
window = freq_window*nchan
window = int(np.ceil(window) // 2*2 +1)

### PEAK OF FITTED POLYNOMIALS / CENTRE OF PULSE BANDWIDTHS
freq_i = []
# centre of bandwidth of sub-pulse
bw_centre = []
for i in range(len(peaks)):
	fit_spec = savgol_filter(eval("S"+str(i)), window, freq_poly)

	# emission bandwidth
	peak_idx = np.argmax(fit_spec)
	dat_freq = np.arange(704.5,4032.5,1)
	spec_peak = dat_freq[peak_idx]
	spec_idx1, spec_idx2, freq_width = cal_fwtm (dat_freq, fit_spec)
	# centre bandwidth
	bw_c = (np.arange(nchan)[spec_idx1]+np.arange(nchan)[spec_idx2])/2
	bw_centre.append(bw_c)

	# index peak of fitted polynomial
	idx = np.argmax(fit_spec)
	freq_i.append(idx)

freq_i = np.array(freq_i)
bw_centre = np.array(bw_centre)


### DYNAMIC SPECTRA ###
if sys.argv[2] == "I":
	c1 = a.clone()
	c1.remove_baseline()
	c1.tscrunch()
	c1.pscrunch()
	data2 = c1.get_data()
	nsub, npol, nchan, nbin = data2.shape

	# on-pulse phase bin start and finish
	ps = int(np.round(p1*nbin))
	pf = int(np.round(p2*nbin))

	# intensity
	I = data2[0,0,:,ps:pf]

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
minorLocator = MultipleLocator(nchan/n_ytick/7.)
ax2.yaxis.set_minor_locator(minorLocator)
ax2.set_yticks(np.linspace(0.0, nchan, n_ytick))
freq0 = a.get_centre_frequency() - a.get_bandwidth()/2.
cfreq = a.get_bandwidth()/nchan
ax2.set_yticklabels([])
ax2.tick_params(axis="y", which='both', direction="in", pad=-22)

ax2.plot(spectrum, np.arange(nchan), ls='-', color='k', lw=2)
ax2.plot(x, np.arange(nchan), ls='--', color='k')


### PLOT DYNAMIC SPECTRUM
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
#fluxes_i, fluxes = zip(*sorted(zip(fluxes_i, fluxes)))

ax1.imshow(eval(p), cmap="Spectral", vmin=np.min(eval(p)), vmax=0.3*np.max(eval(p)), aspect='auto', origin='lower')
# plot peaks on dynamic spectra
#ax1.scatter(fluxes_i, freq_i, marker='x', c='k')
ax1.scatter(peaks, bw_centre, marker='.', c='k')
# plot best fit line
drift_rate_line = np.poly1d(np.polyfit(peaks, bw_centre, 1))(np.unique(peaks))
slope, intercept = np.polyfit(peaks, bw_centre, 1)
ax1.plot(np.unique(peaks), drift_rate_line, linestyle='-', c='k', label=r'$\frac{d\nu}{dt}$ = %s MHz s$^{-1}$'%(np.round(slope,2)))
ax1.legend(loc='lower left')

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
ax0.plot(np.arange(nbin), 0*np.arange(nbin), ls='--', color='k')
# plot peaks and mins
ax0.plot(peaks, flux[peaks], 'x', c='b', label='Peaks > 3 Jy')
ax0.plot(peak_mins, flux[peak_mins], 'x', c='r')

ax0.tick_params(axis="x", which='both', direction="in", pad=-22)
ax0.set_xlim(0.0, pf-ps-1)
ax0.set_xticks(xticks_x)
ax0.set_xticklabels([])

ax0.set_ylim(-5,peak_flux/1000+5)
ax0.set_ylabel('Flux density (Jy)', fontsize=12)

ax0.legend()
ax0.set_title('%s Polarisation %s'%(p,sys.argv[1].split('.')[0]), fontsize=12)


### SAVE FIGURE
plt.savefig("pulse_drift_%s.pdf"%(sys.argv[1].split('.')[0]))
print("pulse_drift_%s.pdf"%(sys.argv[1].split('.')[0]))



