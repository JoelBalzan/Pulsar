import psrchive
import os
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as cm

#########################
### GIANT MICROPULSES ###
#########################

#project code
PCODE = sys.argv[1]
files = sorted(glob.glob("*.rescaled"))
if os.path.isfile('pulse_profiles_'+PCODE+'.npy'):
    pulse_profiles = np.load('pulse_profiles_'+PCODE+'.npy')
else:
    pulse_profiles = []
    counter = 0
    for ar in files:
        a = psrchive.Archive_load(ar)
        a.remove_baseline()
        a.tscrunch()
        a.fscrunch()
        a.pscrunch()
        a.centre()
        data = a.get_data()
        pulse_profiles.append(data[0,0,0,:]/1000)

        # file progress counter
        counter += 1
        print("%s/%s"%(counter,len(files))," files completed", end='\r')
    ## SAVE FLUXES TO FILE
    np.save('pulse_profiles_'+PCODE+'.npy', pulse_profiles)
nfile, nbin = np.shape(pulse_profiles)


### PLOTTING ###
# load first file to get period
a = psrchive.Archive_load(files[0])
# s per bin
period = a.integration_length()
spb = period/nbin

# mean of pulse profiles
mean_profile = np.mean(pulse_profiles, axis=0)

# find brightest bin in each profile
brightest = np.amax(pulse_profiles, axis=0)
# brightest profile bins/mean profile
brightest_mean = brightest/np.max(mean_profile)

# phase start and finish
ps = 0#int(np.argmax(mean_profile)-0.12*nbin)
pf = 2**15#int(np.argmax(mean_profile)+0.12*nbin)
# xticks
xtickslabels = np.round(np.linspace(ps*spb, pf*spb, 11),3)
xticks = np.linspace(0, pf-ps-1, len(xtickslabels))

# trim arrays
brightest_mean = brightest_mean[ps:pf]
mean_profile = mean_profile[ps:pf]
nbin = np.shape(brightest_mean)[0]

A4x, A4y = 8.27, 11.69
# so we can use \mathrm{} in labels
mpl.rcParams['text.usetex'] = False
gmps = []
gmps_i = []
for i in range(nbin):
    if brightest_mean[i] > 50:
        gmps.append(brightest_mean[i])
        gmps_i.append(i)
if len(gmps) != 0:
    #colours = cm.tab20(np.linspace(0, 1, len(gmps)))
    fig = plt.figure(figsize=(A4y,A4x/2), dpi=300)
    scale = np.max(brightest_mean)/np.max(mean_profile)
    plt.plot(mean_profile*scale, c='k', lw=0.5)
    plt.plot(brightest_mean, c='b', lw=0.5, linestyle='--')
    plt.margins(x=0)
    plt.xticks(xticks, xtickslabels, fontsize=10)
    plt.xlabel('Time (s)', fontsize=10)
    plt.ylabel(r'E$_{{i}}/\langle{{\mathrm{E}_{{p}}}}\rangle$', fontsize=10)
    plt.plot(gmps_i[np.argmax(gmps)], np.max(gmps), marker='*', c='r', label=r'${:.1f}\langle\mathrm{{E}}_{{{:.0f}}}\rangle$'.format(np.max(gmps),gmps_i[np.argmax(gmps)]+ps))
    plt.legend()
    plt.savefig('GMP_'+PCODE+'.pdf', dpi=300, bbox_inches='tight')
    print('GMP_'+PCODE+'.pdf')
else:
    print("No GMPs found")

fig = plt.figure(figsize=(A4y,A4x), dpi=300)
xtickslabels = np.round(np.linspace(ps*spb, pf*spb, 11),2)
xticks = np.linspace(0, pf-ps-1, len(xtickslabels))
plt.xlabel('Time (s)')
plt.ylabel('Flux Density (Jy)')
plt.plot(mean_profile, c='k', lw=0.5)
plt.xticks(xticks, xtickslabels)
plt.margins(x=0)
plt.savefig('mean_profile_%s.pdf'%PCODE, dpi=300, bbox_inches='tight')