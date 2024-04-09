import matplotlib
matplotlib.rcParams["font.size"] = 15
matplotlib.rcParams["xtick.labelsize"] = 12 
matplotlib.rcParams["ytick.labelsize"] = 12 

from copy import deepcopy
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt
#import fitburst as fb
from fitburst.analysis.model import SpectrumModeler
import numpy as np
import sys
import psrchive
import glob
import os


a = psrchive.Archive_load(sys.argv[1])

a.remove_baseline()
a.tscrunch()
a.pscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

pulse_profile = np.mean(data[0,0,:,:], axis=0)/1000

bw = a.get_bandwidth()
cf = a.get_centre_frequency()
max_freq = cf+bw/2
min_freq = cf-bw/2

# neutron star period (s)
period = a.integration_length()

# minimum height of peaks (Jy)
h=3
peaks, _ = find_peaks(pulse_profile, height=h, distance=8)
widths = peak_widths(pulse_profile, peaks, rel_height=0.8)


freqs = np.arange(min_freq, max_freq+1, bw/nchan)
times = np.arange(0,period+period/nbin, period/nbin)

# adjust DM value to zero offset, if necessary.
num_components = len(peaks)

# indicate whether the spectrum is de-dispersed or not.
is_dedispersed = True

# define dictiomary containing parameter values.
burst_parameters = {                                                     
    "amplitude"            : list( np.log10(pulse_profile[peaks]) ),    # a list containing the the log (base 10) of the overall signal amplitude
    "arrival_time"         : list(peaks*(period/nbin)),                 # a list containing the arrival times, in seconds
    "burst_width"          : list(widths[0]),                           # a list containing the temporal widths, in seconds
    "dm"                   : [178.5] * num_components,                  # a list containing the dispersion measures (DM), in parsec per cubic centimeter
    "dm_index"             : [-2.] * num_components,                    # a list containing the exponents of frequency dependence in DM delay
    "ref_freq"             : [4032.] * num_components,                  # a list containing the reference frequencies for arrival-time and power-law parameter estimates, in MHz (held fixed)
    "scattering_index"     : [-4.] * num_components,                    # a list containing the exponents of frequency dependence in scatter-broadening
    "scattering_timescale" : [0.] * num_components,                     # a list containing the scattering timescales, in seconds
    "spectral_index"       : [0.] * num_components,                     # a list containing the power-law spectral indices
    "spectral_running"     : [-100.] * num_components,                  # a list containing the power-law spectral running
}


if is_dedispersed:
    burst_parameters["dm"] = [0.] * num_components

# now instantiate the SpectrumModeler for a n-component model.
model = SpectrumModeler(freqs, times, num_components = num_components)

# now update Gaussian-SED model object to use the above values.
model.update_parameters(burst_parameters)

if os.path.isfile("spectrum_model_%s.npy"%(sys.argv[1].split(os.extsep, 1)[0])):
    spectrum_model = np.load("spectrum_model_%s.npy"%(sys.argv[1].split(os.extsep, 1)[0]))
else:
    # grab the model spectrum.
    spectrum_model = model.compute_model()
    np.save("spectrum_model_%s.npy"%(sys.argv[1].split(os.extsep, 1)[0]), spectrum_model)
    print("spectrum_model_%s.npy"%(sys.argv[1].split(os.extsep, 1)[0]))

plt.imshow(spectrum_model, aspect="auto", origin="lower", cmap="viridis")
plt.savefig("spectrum_model_%s.pdf"%(sys.argv[1].split(os.extsep, 1)[0]))
print("spectrum_model_%s.pdf"%(sys.argv[1].split(os.extsep, 1)[0]))