from astropy.io import fits
import numpy as np
import sys
import psrchive

a = psrchive.Archive_load(sys.argv[1])

wts = a.get_weights()
nchan = wts.shape[1]

zapped = []
for i in range(nchan):
    if (wts[0][i] == 0.):
        zapped.append(i)


print('paz -z "%s" -e %s.pazi %s'%((" ".join(str(x) for x in zapped)), sys.argv[1].split('.')[-1], sys.argv[1]))
