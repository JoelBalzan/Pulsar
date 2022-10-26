from astropy.io import fits
import sys


hdul = fits.open(sys.argv[1])
data = hdul[2].data

zapped = []
for i in range(data['DAT_WTS'].shape[1]):
    if (data['DAT_WTS'][0,i] == 0):
        zapped.append(i)


print('paz -z "%s" -e pazi %s'%((" ".join(str(x) for x in zapped)),sys.argv[1]))
