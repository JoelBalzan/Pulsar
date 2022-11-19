import numpy as np
import sys
import subprocess
import glob
import os
import matplotlib.pyplot as plt
import itertools

#all_files = glob.glob("*.rescaled")
#for filename in all_files:
#    cmd = 'pdv -Kx -B 512 %s > pol_frac.txt'%filename

cmd = 'pdv -x -B 512 %s > pol_frac.txt && cut -f2- pol_frac.txt > pf.txt && rm -f pol_frac.txt'%sys.argv[1]
print(cmd)
subprocess.call(cmd, shell=True)

#with open('pol_frac.txt') as f_in:
#    pol = np.genfromtxt(itertools.islice(f_in,1,None,3),dtype=int)

pol = np.genfromtxt('pf.txt',skip_header=1)


x = np.arange(len(pol))

plt.figure(figsize=(15,10),dpi=300)
plt.plot(x, pol[:,4], label='L/I')
plt.plot(x, pol[:,5], label='V/I')
plt.plot(x, pol[:,6], label='|V|/I')
plt.xlabel('Phase bin')
plt.ylabel('Polarisation Fraction')
plt.legend()
plt.savefig("pf_test.png")


