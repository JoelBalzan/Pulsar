import numpy as np
import sys
import psrchive

a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.pscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape

# Phase zoom factor
z = 0.0001

# peak and index
peak_idx = np.argmax(np.mean(data[0,0,:,:], axis=0))
# on-pulse phase start and finish
p1 = np.round(peak_idx/nbin - z, 4)
p2 = np.round(peak_idx/nbin + z, 4)
# on-pulse phase bin start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))


## low-freq peak
#low_freq_peak = np.argmax(np.mean(data[0,0,93:121,:], axis=0))
#
## high-freq peak
#high_freq_peak = np.argmax(np.mean(data[0,0,142:,:], axis=0))
#
## milliseconds per bin
#period = a.integration_length()
#bs = 1000*period/nbin
#
#print(f'Subband offset = {np.round(bs*(high_freq_peak-low_freq_peak), 4)} ms')


# subband 1 peak
subband1_peak = np.argmax(np.mean(data[0,0,0:128,ps:pf], axis=0))

# subband 2 peak
subband2_peak = np.argmax(np.mean(data[0,0,129:257,ps:pf], axis=0))

# peak of rest of the subbands
high_freq_peak = np.argmax(np.mean(data[0,0,258:,ps:pf], axis=0))

# milliseconds per bin
period = a.integration_length()
bs = 1000*period/nbin

print(f'Subband 1 offset = {np.round(bs*(high_freq_peak-subband1_peak), 4)} ms')
print(f'Subband 2 offset = {np.round(bs*(high_freq_peak-subband2_peak), 4)} ms')