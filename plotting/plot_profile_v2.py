import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import gridspec
import argparse
from scipy.signal import savgol_filter
import psrchive

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

#################################################################
####### Main ########
parser = argparse.ArgumentParser(description='Create PRESTO processing database')
# required
parser.add_argument('-i',        '--input_filename',      metavar='filename',      required=True, help='Input file name')
parser.add_argument('-f',        '--avefreq',             metavar='2.0',           required=True, help='Frequency scrunch by a factor')
parser.add_argument('-t',        '--tsamp',               metavar='64',            required=True, help='Time resolution of input file in microsecond (not this is likely to be different from the header)')
# not required
parser.add_argument('-b',        '--avephase',            metavar='1.0',           help='Time scrunch by a factor', default = 1, type = int)
parser.add_argument('-off',      '--offpulse',            metavar='0.2',           help='Off-pulse fraction', default = 0.5, type = float)
parser.add_argument('-o',        '--output_filename',     metavar='filename',      help='Output file name', default = '', type = str)
#parser.add_argument('-ncpus',    '--num_cpus',       metavar='10',             help='Number of CPU to use', default = 10, type = int)
#parser.add_argument('-skipRFI',  action='store_true', default=False,           help='Skip rfifind?')
#parser.add_argument('-skipFold',  action='store_true', default=False,           help='Skip prepfold?')
parser.add_argument('-bw',        '--phase_window',       metavar='0.3',           help='Phase window size', default = 0.3, type = float)
parser.add_argument('-fw',        '--freq_window',       metavar='0.3',           help='Frequency window size', default = 0.3, type = float)
parser.add_argument('-bp',        '--phase_poly',       metavar='5',           help='Phase window poly order', default = 5, type = float)
parser.add_argument('-fp',        '--freq_poly',       metavar='5',           help='Frequency window poly order', default = 5, type = float)
parser.add_argument('-br',        '--phase_range',     metavar='0.3 0.6',  nargs='+',    help='Fix phase window range', default = [0., 0.], type = float)
parser.add_argument('-fr',        '--freq_range',     metavar='100 200',  nargs='+',    help='Fix frequency window range (in channel number)', default = [0., 0.], type = float)
args = parser.parse_args()

filename = args.input_filename
print (filename)
if args.output_filename == '':
	outname = filename.split('.')[0] + '.png'
	result = filename.split('.')[0] + '.txt'
else:
	outname = args.output_filename + '.png'
	result = args.output_filename + '.txt'

ave_freq = int(args.avefreq)
ave_time = int(args.avephase)
tsamp = float(args.tsamp)*1e-6   # second
off = float(args.offpulse)
phase_window = float(args.phase_window)
freq_window = float(args.freq_window)
phase_poly = int(args.phase_poly)
freq_poly = int(args.freq_poly)

phase_idx1 = float(args.phase_range[0])
phase_idx2 = float(args.phase_range[1])
freq_idx1 = int(args.freq_range[0])
freq_idx2 = int(args.freq_range[1])

###########################################################

#filename = sys.argv[1]  # file name of the fits file
#fig = plt.figure(figsize=(10,8))
#ax = fig.add_subplot(211)

fig = plt.figure(figsize=(10, 12), dpi=300) 
gs = gridspec.GridSpec(5, 5, hspace=0, wspace=0, left=0.17) 
ax0 = plt.subplot(gs[0,:4])
ax1 = plt.subplot(gs[1:5, :4])
ax2 = plt.subplot(gs[1:5, 4])
ax0.tick_params(axis='both', which='major', labelsize=12)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)

###########################################################

ar = psrchive.Archive_load(filename)
#print (tsamp, ar.tsamp, ar.tsamp*4)
#print ('WTS', ar.dat_wts.shape, ar.dat_wts)
#ar.get_profile()
ar.remove_baseline()
#ar.remove_baseline(sub=0, delta=off)
ar.tscrunch()
ar.pscrunch()
ar.fscrunch(ave_freq)
if ave_time != 1:
	ar.bscrunch(ave_time)
	tsamp *= ave_time

data = ar.get_data()
nsub, npol, nchan, nbin = data.shape
profile = data[0,0,:,:]
labeltxt = "{0:0.2f}".format(ar.get_mjds()[0])

print ('Nbin: %d; Nchn: %d'%(nbin, nchan))
###########################################################
n_xtick = 6
minorLocator = MultipleLocator(nbin/n_xtick/5.)
ax1.xaxis.set_minor_locator(minorLocator)
ax1.set_xlim(0.0, nbin-1)
ax1.set_xticks(np.linspace(0.0, nbin, n_xtick))
ax1.set_xticklabels(np.round(np.linspace(0, nbin, n_xtick)*tsamp*1000., 1), fontsize=12)
ax1.set_xlabel('Time (ms)', fontsize=14)

n_ytick = 8
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
minorLocator = MultipleLocator(nchan/n_ytick/7.)
ax1.yaxis.set_minor_locator(minorLocator)
ax1.set_ylabel('Frequency (MHz)', fontsize=12)
ax1.set_ylim(0, nchan-1)
ax1.set_yticks(np.linspace(0.0, nchan, n_ytick))
freq0 = ar.get_centre_frequency() - ar.get_bandwidth()/2.
cfreq = ar.get_bandwidth()/nchan
ax1.set_yticklabels(np.round(freq0 + np.linspace(0.0, nchan, n_ytick)*cfreq, 0), fontsize=12)

ax1.imshow(profile, label=labeltxt, aspect='auto', origin='lower', cmap='Spectral', interpolation='none', vmin=np.min(profile), vmax=np.max(profile))
#ax[i,0].plot(phase, vm, ls='-', color='green')
#ax[i,0].legend(loc='upper right', fontsize=5)
#ax[i,0].text(0.65,40,labeltxt,fontsize=9)

print(np.amin(profile), np.amax(profile))

################
ax0.set_title(filename.split('.')[0], fontsize=12)

print ('nchan, nbin: ', profile.shape)
if freq_idx1 == 0. and freq_idx2 == 0.:
	profile_1d = np.mean(profile, axis=0)
else:
	print('Using freqency from {0} to {1}'.format(1600+ave_freq*freq_idx1, 1600+ave_freq*freq_idx2))
	profile_1d = np.mean(profile[freq_idx1:freq_idx2,:], axis=0)

ax0.set_xlim(0.0, nbin-1)
ax0.tick_params(axis="x", which='both', direction="in", pad=-22)
minorLocator = MultipleLocator(nbin/n_xtick/5.)
ax0.xaxis.set_minor_locator(minorLocator)
ax0.set_xticks(np.linspace(0.0, nbin, n_xtick))
ax0.set_xticklabels([])
#ax0.set_yticks([])
phase = np.linspace(0,1,nbin)
y = 0*np.arange(nbin)

ax0.plot(np.arange(nbin), y, ls='--', color='k')
ax0.set_ylabel('Flux density (mJy)', fontsize=12)

if phase_idx1 == 0. and phase_idx2 == 0.:
	window = phase_window*nbin
	window = int(np.ceil(window) // 2*2 +1)
	fit_prof = savgol_filter(profile_1d, window, phase_poly)
	phase_idx1, phase_idx2, prof_width = cal_fwtm (np.arange(nbin), fit_prof)
	print ('phase range: ', phase_idx1, phase_idx2)
else:
	phase_idx1 = int(phase_idx1*nbin)
	phase_idx2 = int(phase_idx2*nbin)
	phase_poly = 3
	phase_window = float((phase_idx2-phase_idx1)/nbin)

#ax0.vlines(np.arange(nbin)[phase_idx1], np.amin(profile_1d), np.amax(profile_1d), ls='--', color='b')
#ax0.vlines(np.arange(nbin)[phase_idx2], np.amin(profile_1d), np.amax(profile_1d), ls='--', color='b')
################
#print (profile.shape, nbin, apeak)
spectrum = np.mean(profile[:,phase_idx1:phase_idx2], axis=1)
x = 0*np.arange(nchan)

ax2.set_ylim(0.0, nchan-1)
ax2.set_yticks(np.arange(0, nchan, 16))
ax2.set_yticklabels([])

n_ytick = 8
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
minorLocator = MultipleLocator(nchan/n_ytick/7.)
ax2.yaxis.set_minor_locator(minorLocator)
ax2.set_yticks(np.linspace(0.0, nchan, n_ytick))
freq0 = ar.get_centre_frequency() - ar.get_bandwidth()/2.
cfreq = ar.get_bandwidth()/nchan
ax2.set_yticklabels([])
ax2.tick_params(axis="y", which='both', direction="in", pad=-22)

ax2.plot(spectrum, np.arange(nchan), ls='-', color='k', lw=2)
ax2.plot(x, np.arange(nchan), ls='--', color='k')

window = freq_window*nchan
window = int(np.ceil(window) // 2*2 +1)
#mask = spectrum != 0.0
#fit_spec = savgol_filter(spectrum[mask], window, freq_poly)
#ax2.plot(fit_spec, np.arange(nchan)[mask], ls='--', color='r')
fit_spec = savgol_filter(spectrum, window, freq_poly)
ax2.plot(fit_spec, np.arange(nchan), ls='--', color='r')

peak_idx = np.argmax(fit_spec)
#spec_peak = ar.dat_freq[0, mask][peak_idx]
dat_freq = np.array([ar.get_first_Integration().get_centre_frequency(i) for i in range(nchan)])
spec_peak = dat_freq[peak_idx]

#spec_idx1, spec_idx2, freq_width = cal_fwtm (ar.dat_freq[0, mask], fit_spec)
spec_idx1, spec_idx2, freq_width = cal_fwtm (dat_freq, fit_spec)
if freq_width != 0.:
	#ax2.hlines(np.arange(nchan)[mask][spec_idx1], np.amin(spectrum), np.amax(spectrum), ls='--', color='r')
	#ax2.hlines(np.arange(nchan)[mask][spec_idx2], np.amin(spectrum), np.amax(spectrum), ls='--', color='r')
	ax2.hlines(np.arange(nchan)[spec_idx1], np.amin(spectrum), np.amax(spectrum), ls='--', color='r')
	ax2.hlines(np.arange(nchan)[spec_idx2], np.amin(spectrum), np.amax(spectrum), ls='--', color='r')

####################
profile_1d = np.mean(profile[spec_idx1:spec_idx2,:], axis=0)
ax0.plot(np.arange(nbin), profile_1d, ls='-', color='k', lw=2)

window = phase_window*nbin
window = int(np.ceil(window) // 2*2 +1)
fit_prof = savgol_filter(profile_1d, window, phase_poly)
phase_idx1, phase_idx2, prof_width = cal_fwtm (np.arange(nbin), fit_prof)
ax0.plot(np.arange(nbin), fit_prof, ls='--', color='r')

ax0.vlines(np.arange(nbin)[phase_idx1], np.amin(profile_1d), np.amax(profile_1d), ls='--', color='r')
ax0.vlines(np.arange(nbin)[phase_idx2], np.amin(profile_1d), np.amax(profile_1d), ls='--', color='r')

peak = np.amax(profile_1d)
apeak = np.argmax(profile_1d)
s = np.sum(profile_1d[phase_idx1-10:phase_idx2+10])
w = s/peak

xb = np.arange(apeak-w/2., apeak+w/2., 0.1)
yb = peak + 0*xb
ax0.plot(xb, yb, ls='--', color='red')

plt.savefig(outname)
#plt.show()
####################

#mjd = ar.stt_imjd + (ar.stt_smjd + ar.stt_offs)/86400.
#
#f = open(result, 'w+')
##f.write('Peak (mJy)    Width (ms)    Fluence (mJy ms)   Freq (MHz)   Width_freq (MHz) MJD\n')
#f.write('{0:.2f}   {1:.2f}    {2:.2f}  {3:.2f}  {4:.2f} {5:.10f} {6}\n'.format(peak, w*tsamp*1000., s*tsamp*1000., spec_peak, freq_width, mjd, filename))
#f.close()
