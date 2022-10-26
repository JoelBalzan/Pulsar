import glob
import numpy as np
import subprocess
import os
from multiprocessing import Pool
import argparse
from astropy.io import fits
import re

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
#parser.add_argument('-ds',      '--dm_step',    metavar='0.1',      required=True, help='DM step')
parser.add_argument('-ncpu',    '--num_cpu',    metavar='5',        required=True, help='Number of CPUs')
parser.add_argument('-e',       '--ext',        metavar='0to3',     required=True, help='Input file extension')
# not required
parser.add_argument('-r',  '--resolution',      metavar='0.512',    help='Acceptable time resolution (ms)', default=0.512, type = float)
parser.add_argument('-cDM',  '--dm_co',         metavar='0.0',      help='Coherent de-dispersion DM', default=0.0, type = float)
parser.add_argument('-nsubband', '--num_subband_dm',      metavar='64',    help='Number of subbands for prepsubband', default=64, type = int)
##########################
args = parser.parse_args()
dmlo = float(args.dm_low)
dmhi = float(args.dm_high)
#dms = float(args.dm_step)
ncpus = int(args.num_cpu)
nsubDM = args.num_subband_dm        # number of subbands for prepsubband
ext = args.ext
resolution = args.resolution
cDM = args.dm_co

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
                cmd = 'python /DATA/CARINA_6/dai02a/UWL/search/single_pulse/run_rfifind_multiFiles.py -f %s -t_rfi 2.0 -zapFreq /DATA/CARINA_6/dai02a/UWL/search/single_pulse/uwl_zap_list -freqsig 4 -timesig 10 -out rfi_root'%(filename)
                subprocess.call(cmd, shell=True)

            ############################
            ########### read in header information ###############
            print ('\n')
            print ('Reading header information from %s\n'%filename)
            hdu = fits.open(filename)
            obsbw = hdu[0].header['OBSBW']             # Total bandwidth (MHz)
            obsfreq = hdu[0].header['OBSFREQ']         # Central frequency (MHz)
            nchn = hdu['SUBINT'].header['NCHAN']       # Total number of channels
            npol = hdu['SUBINT'].header['NPOL']        # Number of polarisation

            nsub = hdu['SUBINT'].header['NAXIS2']      # Number of subint
            nsblk = hdu['SUBINT'].header['NSBLK']      # Samples/row
            nout = nsub*nsblk                # Spectra per file
            tsamp = hdu['SUBINT'].header['TBIN']       # Sampling timie (s)

            hdu.close()

            print ('Total Bandwidth: {0} MHz'.format(obsbw))
            print ('Central frequency: {0} MHz'.format(obsfreq))
            print ('Number of channels: {0}'.format(nchn))
            print ('Sample time: {0} s'.format(tsamp))
            print ('Spectra per file: {0}'.format(nout))

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
                    # Step 2: plan operations for prepsubbands, fft, rednoise, search...
                    ddplan_fig = 'ddplan_fig.eps'
                    cmd = 'DDplan.py -l {0} -d {1} -f {2} -b {3} -n {4} -t {5} -s {6} -o {7} -c {8} -r {9}'.format(dmlo, dmhi, obsfreq, obsbw, nchn, tsamp, nsubDM, ddplan_fig, cDM, resolution)
                    print (cmd)
                    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                    lines = proc.stdout.readlines()

                    dmlos = []
                    ddm = []
                    downsamp = []
                    dm_call = []
                    call = []
                    for line in lines:
                        numbers = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", line.decode('utf-8'))
                        if len(numbers) == 9:
                            dmlos.append(float(line.split()[0]))
                            ddm.append(float(line.split()[2]))      # DM step
                            downsamp.append(int(line.split()[3]))   # down sample
                            dm_call.append(int(line.split()[6]))    # DMs/call
                            call.append(int(line.split()[7]))       # Total number of call

                    cmd_list = []
                    for i in range(len(dmlos)):
                        lodm = dmlos[i]
                        dmstep = ddm[i]
                        for j in range(call[i]):
                            lodm_temp = lodm + j*dm_call[i]*dmstep
                            #dedispersion = '{0} {1} {2} {3} {4} {5} {6}'.format('Prepsubband', lodm_temp, dmstep, dm_call[i], self.nout, self.nsubDM, downsamp[i])
                            cmd_list.append('prepsubband -psrfits -noscales -nooffsets -noweights -ncpus 1 -mask rfi_root_rfifind.mask -lodm {0} -dmstep {1} -numdms {2} -downsamp {3} -numout {4} -nsub {5} -o {6} {7}'.format(lodm_temp, dmstep, dm_call[i], downsamp[i], int(nout/downsamp[i]), nsubDM, filename, filename))

                    #  Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract
                    #    200.000   3008.000    0.50       8   12.00    5616      24     234    1

                    '''
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
                    '''

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
