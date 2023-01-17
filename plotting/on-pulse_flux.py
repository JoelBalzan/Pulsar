import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob
import os


a = psrchive.Archive_load(sys.argv[1])

# load data
c1 = a.clone()
c1.remove_baseline()
c1.tscrunch()
c1.fscrunch()
c1.pscrunch()

data1 = c1.get_data()
nsub, npol, nchan, nbin = data1.shape

# on-pulse start and finish phase bins
on_s = int(float(sys.argv[2])*nbin)
on_f = int(float(sys.argv[3])*nbin)

# on-pulse mean flux density (Jy)
flux = np.sum(data1[0,0,0,on_s:on_f]/1000)/(on_f-on_s)
print("Average flux density = ", np.round(flux,3), " Jy")
