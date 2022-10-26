import glob
import numpy as np
import subprocess
import os
from multiprocessing import Pool
import argparse

def run_job(cmd):
    print (cmd)
    subprocess.call(cmd, shell=True)

#################################################################
####### Main ########
parser = argparse.ArgumentParser(description='Create PRESTO processing database')
# required
parser.add_argument('-s0',      '--band_low',     metavar='0',       required=True, help='Low band')
parser.add_argument('-s1',      '--band_high',    metavar='3',      required=True, help='High band')
parser.add_argument('-ncpu',    '--num_cpu',    metavar='5',        required=True, help='Number of CPUs')
# not required
#parser.add_argument('-zapFreq',  '--zapFreq_list',      metavar='zapFreqFile',    help='File containing all the frequencies to be zapped (MHz)', default = '', type = str)
##########################
args = parser.parse_args()
low = int(args.band_low)
high = int(args.band_high)
ncpus = int(args.num_cpu)

##########################

all_sf = glob.glob('*.sf')
ext = '.%d'%low + 'to' + '%d'%high

cmd_list = []
for filename in all_sf:
    outfile = filename + ext
    if os.path.isfile(outfile):
        print ('%s done!'%outfile)
    else:
        cmd = 'python /DATA/CARINA_6/dai02a/UWL/search/split_search/splitSrch_v2.py -f {0} -o {1} -sub0 {2} -sub1 {3}'.format(filename, outfile, low, high)
        cmd_list.append(cmd)
        #subprocess.call(cmd, shell=True)

##########################
ncpus = 5
pool = Pool(ncpus)
pool.map(run_job, cmd_list)
pool.close()
pool.join()
