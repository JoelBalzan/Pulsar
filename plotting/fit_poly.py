import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import psrchive
from matplotlib.ticker import MultipleLocator
from scipy.signal import savgol_filter
from matplotlib import gridspec

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

	print (idx1, idx2)
	freq1 = freq[idx1]
	freq2 = freq[idx2]
	#print (freq1, freq2)
	bw = np.fabs(freq1-freq2)

	#return freq1, freq2, np.fabs(freq1-freq2)
	return int(idx1), int(idx2), bw	

a = psrchive.Archive_load(sys.argv[1])
c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape

profile = data[0,0,:,:]

p1 = float(sys.argv[2])
p2 = float(sys.argv[3])
ps = int(p1*nbin)
pf = int(p2*nbin)


ticks = np.round(np.linspace(p1,p2,num=11),2)
ticks_x = np.linspace(0,pf-ps+1,num=11)

fig = plt.figure(figsize=(15, 10)) 
ax = plt.subplot(111)
ax.set_xlim(0,pf-ps)
plt.plot(data[0,0,0,ps:pf], zorder=1, label='Flux Density', c='black')
plt.xticks(ticks_x, ticks)
plt.xlabel('Phase')
plt.ylabel('Flux Density')
plt.title('%s'%sys.argv[1].split(os.extsep, 1)[0])


c1 = a.clone()
c.tscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape


freq0 = a.get_centre_frequency() - a.get_bandwidth()/2.
cfreq = a.get_bandwidth()/nchan

phase_window = 0.3
phase_poly = 4
window = phase_window*nbin
window = int(np.ceil(window) // 2*2 +1)
fit_prof = savgol_filter(profile, window, phase_poly)
ps, pf, prof_width = cal_fwtm (np.arange(nbin), fit_prof)
ax.plot(fit_prof, ls='--', color='r')

ax.vlines(np.arange(nbin)[ps], np.amin(profile), np.amax(profile), ls='--', color='r')
ax.vlines(np.arange(nbin)[pf], np.amin(profile), np.amax(profile), ls='--', color='r')

peak = np.amax(profile)
apeak = np.argmax(profile)
s = np.sum(profile[ps-10:pf+10])
w = s/peak

xb = np.arange(apeak-w/2., apeak+w/2., 0.1)
yb = peak + 0*xb
ax.plot(xb, yb, ls='--', color='red')

plt.savefig("fit_poly_%s.pdf"%sys.argv[1].split(os.extsep, 1)[0])
print("fit_poly_%s.pdf"%sys.argv[1].split(os.extsep, 1)[0])



#ticks = np.round(np.linspace(p1,p2,num=11),2)
#ticks_x = np.linspace(0,pf-ps+1,num=11)
#
#plt.figure(figsize=(15,10),dpi=300)
#ax = plt.subplot(111)
#ax.set_xlim(0,pf-ps)
#
#
#plt.plot(data[0,0,0,ps:pf], label='Flux Density', c='black')
#plt.xticks(ticks_x, ticks)
#plt.xlabel('Phase')
#plt.ylabel('Flux Density')
#plt.title('%s'%sys.argv[1].split(os.extsep, 1)[0])
#plt.legend()
#
#plt.savefig("flux_pol_%s.pdf"%sys.argv[1].split(os.extsep, 1)[0])
#print("flux_pol_%s.pdf"%sys.argv[1].split(os.extsep, 1)[0])