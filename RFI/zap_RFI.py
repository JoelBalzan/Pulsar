import psrchive
import sys

a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
data = a.get_data()
nsub, npol, nchan, nbin = a.get_data().shape

zap_list = sys.argv[2].split()
for i in range(npol):
    for j in (zap_list):
        for h in range(nsub):
            profile = a.get_Integration(h).get_Profile(i,int(j))
            amps = profile.get_amps()
            amps[:] = 0.


a.unload("{0}".format(str(sys.argv[1])))
print("{0}".format(str(sys.argv[1])))
