import glob
import numpy as np
import subprocess
import os
from multiprocessing import Pool
import argparse

def do_job (cmd):
    print (cmd)
    subprocess.call(cmd, shell=True)

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

##########################

#################################################################
####### Main ########
parser = argparse.ArgumentParser(description='Create PRESTO processing database')
# required
parser.add_argument('-dl',      '--dm_low',     metavar='100',       required=True, help='DM low')
parser.add_argument('-dh',      '--dm_high',    metavar='100',      required=True, help='DM high')
parser.add_argument('-ds',      '--dm_step',    metavar='0.1',      required=True, help='DM step')
parser.add_argument('-ncpu',    '--num_cpu',    metavar='5',        required=True, help='Number of CPUs')
parser.add_argument('-e',       '--ext',        metavar='0to3',     required=True, help='Input file extension')
# not required
#parser.add_argument('-zapFreq',  '--zapFreq_list',      metavar='zapFreqFile',    help='File containing all the frequencies to be zapped (MHz)', default = '', type = str)
##########################
args = parser.parse_args()
dmlow = float(args.dm_low)
dmhigh = float(args.dm_high)
dms = float(args.dm_step)
ncpus = int(args.num_cpu)
ext = args.ext

all_sf = glob.glob('*.%s'%ext)

for filename in all_sf:
    singlepulse = filename + '_singlepulse.ps'
    if os.path.isfile(singlepulse):
        print ('%s done!\n'%singlepulse)
    else:
        dirname = filename.split('.')[0]
        if os.path.exists(dirname):
            print (dirname)
        else:
            os.mkdir(dirname)

        with cd(dirname):
            print (os.getcwd())
            if os.path.isfile(filename):
                print ('%s is here.\n'%(filename))
            else:
                cmd = 'ln -s ../%s .'%(filename)
                subprocess.call(cmd, shell=True)

            ###### run rfi_find #########
            rfi_mask = 'rfi_root_rfifind.mask'
            if os.path.isfile(rfi_mask):
                print ('rfifind done.\n')
            else:
                cmd = 'python /DATA/DRACO_2/dai02a/UWL/search/single_pulse/run_rfifind_multiFiles.py -f %s -t_rfi 2.0 -zapFreq /DATA/DRACO_2/dai02a/UWL/search/single_pulse/uwl_zap_list -freqsig 4 -timesig 10 -out rfi_root'%(filename)
                subprocess.call(cmd, shell=True)

            ###### de-dispersion #########
            dat_list = glob.glob('uwl*.dat')
            inf_list = glob.glob('uwl*.inf')
            pulse_list = glob.glob('uwl*.singlepulse')

            ###### prepsubband #########
            #if len(dat_list) != 0 and len(dat_list) == len(inf_list):
            #    print ('prepsubband done.\n')
            #else:
            #    cmd = 'python ../run_prepsubband_multiCPU.py -f %s -dmlo 1180 -dmhi 1250 -nsubband 64 -ncpus 4 -cDM 1227.7'%(filename)
            #    #cmd = 'python ../run_prepsubband_multiCPU.py -f %s -dmlo 1100 -dmhi 1500 -nsubband 64 -ncpus 4 -cDM 1227.7'%(filename)
            #    print cmd
            #    subprocess.call(cmd, shell=True)

            if (len(pulse_list) == len(inf_list)) and (len(pulse_list) != 0):
                print ('Search done.\n')
            else:
                ##### prepdata #########
                if len(dat_list) != 0 and len(dat_list) == len(inf_list):
                    print ('prepsubband done.\n')
                else:
                    cmd_list = []
                    for dm in np.arange(dmlow, dmhigh, dms):
                    #for dm in np.arange(1743.0,1803.0,0.2):
                    #for dm in np.arange(1190,1280,0.2):
                        inf_name = filename + '_DM%.1f'%dm + '.inf'
                        dat_name = filename + '_DM%.1f'%dm + '.dat'
                        if (inf_name in inf_list) and (dat_name in dat_list):
                            print ('%s done\n'%dat_name)
                        else:
                            cmd_list.append('prepdata -psrfits -noscales -nooffsets -zerodm -nobary -dm {0} -mask rfi_root_rfifind.mask -o {1} {2}'.format(dm, filename + '_DM%.1f'%dm, filename))
                    #print(cmd_list)
                    pool = Pool(ncpus)
                    pool.map(do_job, cmd_list)
                    pool.close()
                    pool.join()

                ###### single pulse search #########
                #cmd = '/DATA/DRACO_2/zha235/ToO_J1708-4009/single_pulse_search.py_new -t 7.0 -b -p -m 1 *.dat'
                all_dat = sorted(glob.glob('*.dat'))
                cmd_list = []
                for datfile in all_dat:
                    #cmd_list.append('/DATA/DRACO_2/zha235/ToO_J1708-4009/single_pulse_search.py_new -t 6.0 -b -p -m 1 %s'%datfile)
                    cmd_list.append('single_pulse_search.py -t 6.0 -b -p -m 1 %s'%datfile)

                pool = Pool(ncpus)
                pool.map(do_job, cmd_list)
                pool.close()
                pool.join()

                #cmd = '/DATA/DRACO_2/zha235/ToO_J1708-4009/single_pulse_search.py_new -t 6.0 -b -p -m 1 *.dat'
                #print cmd
                #subprocess.call(cmd, shell=True)

                #cmd = '/DATA/DRACO_2/zha235/ToO_J1708-4009/single_pulse_search.py_new *.singlepulse'
                cmd = 'single_pulse_search.py *.singlepulse'
                print (cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'mv *_singlepulse.ps ../'
                print (cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'rm -f uwl*.dat'
                print (cmd)
                subprocess.call(cmd, shell=True)
