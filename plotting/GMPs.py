import psrchive
import os
import numpy as np
import glob
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

for i in range(nfile):
    fig = plt.figure(figsize=(10,5), dpi=300)
    gmps = []
    gmps_i = []
    for j in range(nbin):
        if pulse_profiles[i][j]/mean_profile[j] > 100:
            gmps.append(pulse_profiles[i][j]/mean_profile[j])
            gmps_i.append(j)
    if len(gmps) != 0:
        colours = cm.tab20(np.linspace(0, 1, len(gmps)))
        plt.plot(pulse_profiles[i]*np.max(pulse_profiles[i]/mean_profile)/np.max(pulse_profiles[i]), c='k', lw=0.5, label='Pulse Profile')
        plt.plot(pulse_profiles[i]/mean_profile, c='b', lw=0.5, linestyle='--', label='Pulse Profile / Mean Profile')
        plt.xlabel('Phase Bin')
        plt.ylabel('Flux Density (Jy)')
        for k in range(len(gmps)):
            plt.axvline(x=gmps_i[k], c=colours[k], lw=0.5, linestyle='--', label=r'%s$\langle{{E_{{%s}}}}\rangle$'%(gmps[k],gmps_i[k]))
        plt.legend()
        plt.savefig('GMP_'+files[i].split(os.extsep, 1)[0]+'.pdf', dpi=300)
        plt.close()
fig = plt.figure(figsize=(10,5), dpi=300)
plt.plot(mean_profile)
plt.savefig('mean_profile_%s.pdf'%PCODE, dpi=300)