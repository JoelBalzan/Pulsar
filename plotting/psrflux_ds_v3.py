import os
import sys

import psrchive
import numpy as np
import pandas as pd
import glob



PCODE = sys.argv[1]

files = sorted(glob.glob("*calib.rescaled"))
nfiles = len(files)

f = psrchive.Archive_load(files[0])
freq_info = f.get_frequencies()
period = f.integration_length()/60 #minutes
nbins = f.get_nbin()
nchans = f.get_nchan()

MJD = f.get_mjds()[0]

### FREQ ZOOM
bw = f.get_bandwidth()
cf = f.get_centre_frequency()
#lowest observed frequency
min_freq = cf-bw/2
# fscrunching factor
f_scr = bw/f.get_nchan()
del f

isub = np.repeat(np.arange(nfiles), nchans).astype(int)
ichan = np.tile(np.arange(nchans), nfiles).astype(int)
time = np.repeat(np.arange(period/2, nfiles*period, period), nchans)
#freqs = np.tile(np.arange(f1+0.5, f2, f_scr), nfiles)
freqs = np.tile(freq_info, nfiles)
flux = []
flux_err = []
counter = 0
for ar in files:
	a = psrchive.Archive_load(ar)
	a.remove_baseline()
	a.tscrunch()
	a.pscrunch()
	data = a.get_data()
	nsub, npol, nchan, nbin = data.shape

	spectra = np.mean(data[0,0,:,:], axis=1)
	spectra_err = 0.1 * spectra #10% error

	flux.append(spectra)
	flux_err.append(spectra_err)

	# file progress counter
	counter += 1
	print("%s/%s"%(counter,len(files))," files completed", end='\r')

flux = np.array(flux).flatten()
flux_err = np.array(flux_err).flatten()

psrflux = pd.DataFrame({
		'isub': isub,
		'ichan': ichan,
		'time(min)': time,
		'freq(MHz)': freqs,
		'flux': flux,
		'flux_err': flux_err
		})
psrflux.to_pickle("%s_psrflux_v3.pkl"%(PCODE))

np.savetxt("%s_psrflux_v3.ds"%(PCODE), psrflux.values, fmt='%4d %5d %12.4f %14.6f %+14.6e %+14.6e', delimiter='\t', header='\t'.join(psrflux.columns), 
		   comments='# Dynamic spectrum computed by psrflux_ds python script \n# Command line: python %s %s \n# Data file: *.calib.rescaled \n# Flux method: StandardFlux \n# Flux units: Jansky \n# MJD0: %s \n# Data columns: \n# '%(sys.argv[0], sys.argv[1], MJD),
		   )
