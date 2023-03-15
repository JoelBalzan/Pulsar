import numpy as np
import argparse
import matplotlib.pyplot as plt
import psrchive
import os


### ARGUMENTS
parser = argparse.ArgumentParser(description="Generates a Phase vs Frequency plot from an archive file",
								 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-f", "--file", type=str, required=True, help="Archive file")
parser.add_argument("-p", "--polarisation", type=str, required=True, help="Polarisation to plot (I (total intensity after p-scrunching), L (linear), SI, SQ, SU, SV)")
parser.add_argument("-b", "--bscr", type=int, default=0, help="B-scrunch by this factor")
parser.add_argument("-z", "--zoom", type=str, default=(0,1), help="Zoom into this phase range, i.e. 0,0.5")
args = vars(parser.parse_args())

### VARIABLES
pol = args["polarisation"]
b = args["bscr"]
z = str(args["zoom"]).replace('(','').replace(')','')
z = [float(i) for i in z.split(",")]

file = args["file"]
a = psrchive.Archive_load(file)
a.remove_baseline()
a.tscrunch()
if b != 0:
	a.bscrunch(b)

### POLARISATION
if pol == "I":
	a.pscrunch()
	data = a.get_data()
	nsub, npol, nchan, nbin = data.shape

	# phase bin start and finish
	ps = int(np.round(z[0]*nbin))
	pf = int(np.round(z[1]*nbin))

	I = data[0,0,:,ps:pf]
else:
	data = a.get_data()
	nsub, npol, nchan, nbin = data.shape

	# phase bin start and finish
	ps = int(np.round(z[0]*nbin))
	pf = int(np.round(z[1]*nbin))

	# polarisations
	if pol == "SI":
		SI = data[0,0,:,ps:pf]
	if pol == "SQ":
		SQ = data[0,1,:,ps:pf]
	if pol == "SU":
		SU = data[0,2,:,ps:pf]
	if pol == "L":
		L = np.sqrt(data[0,1,:,ps:pf]**2+data[0,2,:,ps:pf]**2)
	if pol == "SV":
		SV = data[0,3,:,ps:pf]


### PLOTTING
plt.figure(figsize=(30,30), dpi=300)

xticks = np.round(np.linspace(z[0],z[1],num=10),4)
xticks_x = np.linspace(0,pf-ps-1,num=len(xticks))
if nchan==208:
	yticks = np.linspace(0,nchan-1, num=nchan).astype(int)
else:
	yticks = np.linspace(0,nchan-1, num=56).astype(int)
yticks_y = np.linspace(0,nchan-1, len(yticks))

ax = plt.subplot(111)
plt.imshow(eval(pol), cmap='viridis', vmin=0.5*np.min(eval(pol)), vmax=0.3*np.max(eval(pol)), aspect='auto', origin='lower', interpolation=None)
plt.xticks(xticks_x, xticks)
plt.yticks(yticks_y, yticks)
plt.xlabel('Phase')
plt.ylabel('Frequency Channel')

### SAVE FIGURE
plt.savefig('phase_freq_%s_%s.pdf'%(pol,file.split(os.extsep, 1)[0]))
print('phase_freq_%s_%s.pdf'%(pol,file.split(os.extsep, 1)[0]))

