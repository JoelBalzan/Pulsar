import numpy as np
import argparse

##########################################################
parser = argparse.ArgumentParser(description='To determine the first sample and total width')
parser.add_argument('-dm', '--cand_dm',  metavar='1200.0',  nargs='+', required=True, help='Candidate DM')
parser.add_argument('-f1', '--freq_low', metavar='704', nargs='+', required=True, help='Frequency lower boundary (MHz)')
parser.add_argument('-f2', '--freq_high', metavar='4032', nargs='+', required=True, help='Frequency higher boundary (MHz)')
parser.add_argument('-s', '--burst',  metavar='120300 1856', nargs='+', required=True, help='Detected burst time (in sample) and highest frequency (MHz)')
parser.add_argument('-dt', '--tsamp', metavar='32.e-6', required=True, help='Sampling time (s)')
######## not required ##########
parser.add_argument('-o',  '--offset', metavar='1000', default=1000, type=int, help='Offset in number of samples')

args = parser.parse_args()
dm = float(args.cand_dm[0])
f1 = float(args.freq_low[0])/1000.
f2 = float(args.freq_high[0])/1000.  # GHz
off = int(args.offset)
tsamp = float(args.tsamp)
s0 = int(args.burst[0])
f0 = float(args.burst[1])/1000.   # GHz

#print (dm, f1, f2, off)
##########################################################
### where is the first sample ########
dt = 4.15*dm*(f0**(-2) - f2**(-2))
samp_start = s0 - dt/1000./tsamp - off

### How many samples in total? ########
dt = 4.15*dm*(f1**(-2) - f2**(-2))
tmp = dt/1000./tsamp + off
width = int((int(tmp/16.) + 5)*16.)

print ('First sample: %d'%samp_start)
print ('Total width: %d'%width)

