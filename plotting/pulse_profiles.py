import psrchive
import os
import numpy as np
import glob

PCODE = "P970"
if os.path.isfile('pulse_profiles_'+PCODE+'.npy'):
    pulse_profiles = np.load('pulse_profiles_'+PCODE+'.npy')
else:
    pulse_profiles = []

    
    counter = 0
    for ar in glob.glob("*.rescaled"):
        a = psrchive.Archive_load(ar)
        a.remove_baseline()
        a.tscrunch()
        a.fscrunch()
        a.pscrunch()
        data = a.get_data()
        pulse_profiles.append(data[0,0,0,:])

        # file progress counter
        counter += 1
        print("%s/%s"%(counter,len(glob.glob("*.rescaled")))," files completed", end='\r')
    ## SAVE FLUXES TO FILE
    np.save('pulse_profiles_'+PCODE+'.npy', pulse_profiles)

