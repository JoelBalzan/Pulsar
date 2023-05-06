import psrchive
import os
import numpy as np
import glob
import matplotlib.pyplot as plt

#########################
### GIANT MICROPULSES ###
#########################

PCODE = "P970"
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

fig = plt.figure(figsize=(10,5), dpi=300)
for i in range(nfile):
    gmps = []
    for j in range(nbin):
        if pulse_profiles[i][j]/mean_profile[j] > 1000:
            gmps.append(j)
    plt.plot(pulse_profiles[i], c='k', lw=0.5, label='Pulse Profile')
    plt.plot(pulse_profiles[i]/mean_profile, c='b', lw=0.5, linestyle='--', label='Pulse Profile / Mean Profile')
    plt.xlabel('Phase Bin')
    plt.ylabel('Flux Density (Jy)')
    for k in range(len(gmps)):
        plt.axvline(x=gmps[k], c='r', lw=0.5, linestyle='--')
    plt.legend()
    plt.savefig('GMP_'+'_'+files[i].split(os.extsep, 1)[0]+'.pdf', dpi=300)
fig = plt.figure(figsize=(10,5), dpi=300)
plt.plot(mean_profile)
plt.savefig('mean_profile.pdf', dpi=300)