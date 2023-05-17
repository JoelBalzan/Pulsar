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
        data = a.get_data()
        pulse_profiles.append(data[0,0,0,:]/1000)

        # file progress counter
        counter += 1
        print("%s/%s"%(counter,len(files))," files completed", end='\r')
    ## SAVE FLUXES TO FILE
    np.save('pulse_profiles_'+PCODE+'.npy', pulse_profiles)
nfile, nbin = np.shape(pulse_profiles)

mean_profile = np.mean(pulse_profiles, axis=0)
ps = int(np.argmax(mean_profile)-0.1*nbin)
pf = int(np.argmax(mean_profile)+0.1*nbin)
xtickslabels = np.round(np.linspace(ps/nbin, pf/nbin, 11),3)
xticks = np.linspace(0, pf-ps-1, len(xtickslabels))

A4x, A4y = 8.27, 11.69
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
max_gmps = []
for i in range(nfile):
    fig = plt.figure(figsize=(A4y,A4x/2), dpi=300)
    gmps = []
    gmps_i = []
    for j in range(nbin):
        if pulse_profiles[i][j]/mean_profile[j] > 100:
            gmps.append(pulse_profiles[i][j]/mean_profile[j])
            gmps_i.append(j)
    if len(gmps) != 0:
        max_gmps.append(np.max(gmps))
        colours = cm.tab20(np.linspace(0, 1, len(gmps)))
        plt.plot(pulse_profiles[i][ps:pf]/mean_profile[ps:pf], c='b', lw=0.5, linestyle='--', label='Pulse Profile/Mean Profile')
        # scale pulse profile to mean profile
        scale = np.max(pulse_profiles[i][ps:pf]/mean_profile[ps:pf])/np.max(pulse_profiles[i][ps:pf])
        plt.plot(pulse_profiles[i][ps:pf]*scale, c='k', lw=0.5, label='Pulse Profile')
        plt.margins(x=0)
        plt.xticks(xticks, xtickslabels)
        plt.xlabel('Phase (Turns)')
        plt.ylabel(r'E$_{{i}}/\langle{{E_{{p}}}}\rangle$')
        for k in range(len(gmps)):
            plt.plot(gmps_i[k], gmps[k]+0.05*gmps[k], c=colours[k], marker='*', label=r'{:.2f}\langle{{E}}_{{:.0f}}\rangle'%(gmps[k],gmps_i[k]))
        plt.legend()
        plt.savefig('GMP_'+files[i].split(os.extsep, 1)[0]+'.pdf', dpi=300, bbox_inches='tight')
    plt.close(fig)
pdfs = sorted(glob.glob("GMP_pulse_*.pdf"))
print("largest GMP:", np.max(max_gmps), "in file:", pdfs[np.argmax(max_gmps)])
fig = plt.figure(figsize=(A4y,A4x), dpi=300)
ps = int(np.argmax(mean_profile)-7000)
pf = int(np.argmax(mean_profile)+3500)
xtickslabels = np.round(np.linspace(ps/nbin, pf/nbin, 11),3)
xticks = np.linspace(0, pf-ps-1, len(xtickslabels))
plt.xlabel('Phase (Turns)')
plt.ylabel('Flux Density (Jy)')
plt.plot(mean_profile[ps:pf], c='k', lw=0.5)
plt.xticks(xticks, xtickslabels)
plt.margins(x=0)
plt.savefig('mean_profile_%s.pdf'%PCODE, dpi=300, bbox_inches='tight')