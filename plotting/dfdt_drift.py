
import numpy as np
import psrchive
import dfdt
import sys
from scipy.signal import find_peaks, peak_widths

a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.pscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

# ms per bin
period = a.integration_length()
mspb = 1000*period/nbin
# folded profile and peak index
F = np.mean(data[0,0,:,:], axis=0)
peak_idx = np.argmax(F)
if sys.argv[1] == 'pulse_65080037.calib.rescaled':
		peaks, _ = find_peaks(F)
		peak_idx = np.where(F==np.sort(F[peaks])[-2])[0][0]
# phase window
width = np.round(20*peak_widths(F, np.array([peak_idx]), rel_height=0.5)[0]).astype(int)
# on-pulse phase bin start and finish
ps = int(np.round(peak_idx - width/2))
pf = int(np.round(peak_idx + width/2))

### FREQ ZOOM
bw = a.get_bandwidth()
cf = a.get_centre_frequency()
#lowest observed frequency
min_freq = cf-bw/2
# fscrunching factor
f_scr = bw/a.get_nchan()
#fs = int((f1-min_freq)/f_scr)
#ff = int((f2-min_freq)/f_scr)
# intensity
dedispersed_intensity = data[0,0,:,ps:pf]

# burst parameters
dm_uncertainty = 0.1  # pc cm-3
source = "XTE 1810-197"
eventid = sys.argv[1]

# instrument parameters
dt_s = period/nbin
df_mhz = a.get_bandwidth()/nchan
bw = a.get_bandwidth()
cf = a.get_centre_frequency()
freq_bottom_mhz = cf-bw/2
freq_top_mhz = bw+freq_bottom_mhz

ds = dfdt.DynamicSpectrum(dt_s, df_mhz, nchan, freq_bottom_mhz, freq_top_mhz)

constrained, dfdt_data, dfdt_mc, dfdt_mc_low, dfdt_mc_high = dfdt.ac_mc_drift(
    dedispersed_intensity, dm_uncertainty, source, eventid, ds,
    dm_trials=100, mc_trials=100
)