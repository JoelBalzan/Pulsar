import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import psrchive
import os
from scipy.signal import find_peaks, peak_widths
from scipy.signal import savgol_filter
from numpy import nan
import math


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

#flux = data[0,0,0,ps:pf]/1000
flux = data[0,0,0,:]/1000

h=3
peaks, _ = find_peaks(flux, height=h)

# peak minimas
mins, _ = find_peaks(-flux)
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
#spectra = data1[0,0,:,ps:pf]/1000
spectra = data1[0,0,:,:]/1000


### DEFINE SPECTRA VARIABLES
dict={}
for i in range(len(peaks)):
	key = str("S"+str(i))
	spec = np.mean(spectra[:,peak_mins[i][0]:peak_mins[i][1]], axis=1)
	dict[key] = spec.tolist()
for key,value in dict.items():
	exec(f'{key}={value}')


#### PLOTTING ####
fig = plt.figure(figsize=(15,15),dpi=300)
g = gridspec.GridSpec(ncols=1, nrows=len(peaks), hspace=0)

freq_window = float(sys.argv[2])
freq_poly = int(sys.argv[3])
window = freq_window*nchan
window = int(np.ceil(window) // 2*2 +1)


for i in range(len(peaks)):
	ax = fig.add_subplot(g[i])
	ax.plot(eval("S"+str(i)), c='k', label='peak %s'%i)
	fit_spec = savgol_filter(eval("S"+str(i)), window, freq_poly)
	ax.plot(np.arange(nchan), fit_spec, ls='--', color='r', label='Deg = %s, Window = %s'%(freq_poly, freq_window))
	ax.plot(np.arange(nchan),0*np.arange(nchan), ls='--', c='k')

	peak_idx = np.argmax(fit_spec)
	dat_freq = np.array([a.get_first_Integration().get_centre_frequency(i) for i in range(nchan)])
	spec_peak = dat_freq[peak_idx]
	spec_idx1, spec_idx2, freq_width = cal_fwtm (dat_freq, fit_spec)
	ax.vlines(np.arange(nchan)[spec_idx1], np.amin(eval("S"+str(i))), np.amax(eval("S"+str(i))), ls='--', color='r')
	ax.vlines(np.arange(nchan)[spec_idx2], np.amin(eval("S"+str(i))), np.amax(eval("S"+str(i))), ls='--', color='r')

	ax.set_xlim(0,nchan-1)
	ax.legend(loc='upper left')
ax.set_xticks(np.linspace(0,nchan-1, 14))
ax.set_xticklabels(np.linspace(704,4032, num=14).astype(int))
ax.set_xlabel('Frequency (MHz)')


### SAVE FIGURE
plt.savefig("sub-pulse_spectra_%s.pdf"%(sys.argv[1].split('.')[0]), bbox_inches='tight')
print("sub-pulse_spectra_%s.pdf"%(sys.argv[1].split('.')[0]))


