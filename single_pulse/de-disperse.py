import math
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from matplotlib import rc
#from matplotlib.patches import Circle

#import subprocess
import argparse
#import glob

#import pyfits
from astropy.io import fits

#rc('text', usetex=True)
# read in parameters

#filename   bw   nchan   freq   tsamp
#uwl_200806_052632_57.sf.wideband   896.000   896   1920.000   3.2e-05

# DM      Sigma      Time (s)     Sample    Downfact
#1215.10   12.37     90.384800    2824525     150

def dedisperse (dat, dm, f1, f2, nchan, tsamp):
    # f1, f2 in GHz
    # tsamp in ms

    chn_bw = (f2 - f1)/nchan

    i = 0
    for f in np.arange(f1 + chn_bw/2., f2, chn_bw):
        dt = 4.15 * dm * (f**(-2.) - f2**(-2.))   # ms
        ds = int(dt/tsamp)    # dispersive delay in number of samples

        #if ds != 0:
        #    print f, f2, ds
        dat[i,:] = np.roll(dat[i,:], -ds) 
        #dat[i,:] = np.roll(dat[i,:], ds) 
        i = i + 1

def zap_rfi (dat, nchan):
    ######### Remove strong narrow band RFI ############
    spec = np.sum(dat, axis=1)
    std = np.std(dat, axis=1)
    n0 = 0
    n1 = 1

    mask = np.ones(nchan, dtype=bool)
    while n0 != n1:
        n0 = nchan - np.count_nonzero(mask)
        std_spec = np.std(spec[mask])
        mean_spec = np.mean(spec[mask])

        std_std = np.std(std[mask])
        mean_std = np.mean(std[mask])
        for i in range(nchan):
            #if np.abs(spec[i]-mean_spec) > 3*std_spec:
            if (np.abs(spec[i]-mean_spec) > 3*std_spec) or (np.abs(std[i]-mean_std) > 3*std_std):
                mask[i] = False

        n1 = nchan - np.count_nonzero(mask)

    dat_mean = np.mean(dat[mask,:])
    m, nbin = dat.shape
    for i in range(nchan):
        if mask[i] == False:
            for j in range(nbin):
                dat[i,j] = dat_mean

    ######### Remove strong wide band RFI ############
    zerodm = np.mean(dat, axis=0)
    for i in range(nchan):
        dat[i,:] = dat[i,:] - zerodm 

def bin_data (dat, num_bin, num_chn):
    m,n = dat.shape
    dat_plot = np.mean(dat.reshape(int(m/num_chn), int(num_chn), n), axis=1)

    m,n = dat_plot.shape
    dat_plot = np.mean(dat_plot.reshape(m, int(n/num_bin), int(num_bin)), axis=2)

    return dat_plot

def cal_intensity (dat, freq, onpulse):
    #### calculate accumulated intensity #########
    dat_use = np.sum(dat[onpulse[0]:onpulse[1],:], axis=0)
    off = np.sum(dat[(onpulse[0]-40):(onpulse[1]-40),:], axis=0)
    num = len(freq)
    intensity = []
    for i in range(num):
        intensity.append(np.mean(dat_use[:i])/float(np.mean(off[:i])))
        #intensity.append(dat_use[i]/float(off[i]))

    return np.array(intensity)

#############################

parser = argparse.ArgumentParser(description='Read PSRFITS format search mode data')
parser.add_argument('-f',  '--input_file',  metavar='Input file name',  nargs='+', required=True, help='Input file name')
parser.add_argument('-dm', '--frb_dm', metavar='DM of the FRB', nargs='+', required=True, help='DM of the FRB')
parser.add_argument('-t',  '--t_sample', metavar='Sampling time', nargs='+', required=True, help='Sampling time (ms)')
parser.add_argument('-f0',  '--freq_start', metavar='Starting frequency (MHz)', nargs='+', required=True, help='Starting frequency')
parser.add_argument('-f1',  '--freq_end', metavar='Ending frequency (MHz)', nargs='+', required=True, help='Ending frequency')
parser.add_argument('-cb',  '--chan_bw', metavar='Channel bandwidth', nargs='+', required=True, help='Channel bandwidth')
parser.add_argument('-fn',  '--num_freq_binning', metavar='Binning', nargs='+', required=True, help='Binning')
parser.add_argument('-bn',  '--num_time_binning', metavar='Binning', nargs='+', required=True, help='Binning')
parser.add_argument('-xlim',  '--xaxis_limit', metavar='Number of sample', nargs='+', required=True, help='Number of sample')
parser.add_argument('-vmin',  '--value_min', metavar='Vmin for imshow', nargs='+', required=True, help='Vmin for imshow')
parser.add_argument('-skip',  action='store_true', default=False,           help='Skip RFI zapping and de-dispersion?')

args = parser.parse_args()
fname = args.input_file[0]
dm = float(args.frb_dm[0])
tsamp = float(args.t_sample[0])
f1 = float(args.freq_start[0])
f2 = float(args.freq_end[0])
chan_bw = float(args.chan_bw[0])
fn = int(args.num_freq_binning[0])
bn = int(args.num_time_binning[0])
xlim = int(args.xaxis_limit[0])
vmin = float(args.value_min[0])
skip = args.skip

nchan = int((f2-f1)/chan_bw)
print (nchan)

###########################################
plt.figure(figsize=(8,10))
dat = np.load(fname)
print (dat.shape)

###########################################
if skip == False:
    dat = np.flip(np.rot90(dat, k=3), axis=1)

    ######### Zap RFI ##########
    print('Zapping...\n')
    zap_rfi (dat, nchan)

    print('De-dispersing...\n')
    dedisperse (dat, dm, f1/1000., f2/1000., nchan, tsamp)
    
    print('Saving...\n')
    tmp = fname.split('.')[0]
    np.save(tmp + '_dedispersed', dat)
else:
    print('Already de-dispersed.\n')

m,n = dat.shape
print (m,n)
prof = np.sum(dat, axis=0)

######################################
gs1 = gridspec.GridSpec(5, 5)
gs1.update(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.01, hspace=0.01)

ax = []
ax.append(plt.subplot(gs1[0:1, 0:4]))
ax.append(plt.subplot(gs1[1:, 0:4]))
ax.append(plt.subplot(gs1[1:, 4]))

######################################
ax[0].set_title('Profile (2020/10/25)')
ax[0].tick_params(axis="x",direction="in", pad=-22)
ax[0].set_xlim(0,int(xlim/float(bn)))
ax[0].set_xticks(np.arange(0,int(xlim/float(bn)),int(xlim/float(bn)/10)))
ax[0].set_xticklabels([])
prof = np.add.reduceat(prof, range(0, n, bn))
n = len(prof)
ymin = np.amin(prof[:int(xlim/float(bn))])
ymax = np.amax(prof[:int(xlim/float(bn))])
ax[0].set_yticks([])
ax[0].set_ylim([ymin, ymax])
ax[0].plot(np.arange(n), prof, color='k')

#ax[1].set_title('De-dispersed')
ax[1].set_xlabel(r'Sample (32$\mu$s)', fontsize=12)
ax[1].set_xlim(0,xlim/float(bn))
ax[1].set_xticks(np.arange(0,int(xlim/float(bn)),int(xlim/float(bn)/10)))
ax[1].set_xticklabels(bn*np.arange(0,int(xlim/float(bn)),int(xlim/float(bn)/10)), fontsize=10)
ax[1].set_ylim([0,nchan/fn])
ax[1].set_yticks(np.arange(0,nchan/fn,nchan/fn/8))
ax[1].set_yticklabels(f1 + np.arange(0,nchan/fn,nchan/fn/8)*fn, fontsize=10)
#ax.imshow(dat, aspect='auto', origin='lower')
dat_plot = bin_data (dat, bn, fn)
print (dat_plot.shape)
print (np.amin(dat_plot), np.amax(dat_plot))
#ax[1].imshow(dat_plot, aspect='auto', origin='lower', vmin=20, cmap='binary')
ax[1].imshow(dat_plot, aspect='auto', origin='lower', vmin=vmin, cmap='binary')

intensity = cal_intensity (dat_plot, np.arange(0,nchan/fn), [120,130])
print (intensity.shape)
ax[2].set_ylim([0,nchan/fn])
ax[2].set_yticks([])
ax[2].set_xticks([])
ax[2].plot(intensity, np.arange(0,nchan/fn), color='k')

plt.show()
