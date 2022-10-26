from astropy.io import fits
import numpy as np
import sys


hdul = fits.open(sys.argv[1])
data = hdul['SUBINT'].data
wts = data['DAT_WTS']

nsub = wts.shape[0]
nchan = wts.shape[1]

zapped = []
for i in range(nchan):
    if ((wts[-1,i] == 0)):
        zapped.append(i)


print('paz -z "%s" -e pazi %s'%((" ".join(str(x) for x in zapped)),sys.argv[1]))
