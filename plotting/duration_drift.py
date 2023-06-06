import glob
import os
import sys

import lmfit
import matplotlib.pyplot as plt
import numpy as np
import psrchive
from lmfit import Model
from lmfit.lineshapes import lorentzian
from lmfit.model import save_modelresult
from matplotlib import gridspec
from matplotlib.patches import Ellipse
from scipy import signal
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_widths


### FUNCTIONS ###
def gauss(x, H, A, x0, sigma):
	return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gauss_fit(x, y):
	mean = sum(x * y) / sum(y)
	sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
	popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
	return popt


def RotGauss2D(x, y, A, x0, y0, sigma_x, sigma_y, theta):
	"""Rotated 2D Gaussian function"""
	theta = np.radians(theta)
	sigx2 = sigma_x**2; sigy2 = sigma_y**2
	a = np.cos(theta)**2/(2*sigx2) + np.sin(theta)**2/(2*sigy2)
	b = np.sin(theta)**2/(2*sigx2) + np.cos(theta)**2/(2*sigy2)
	c = np.sin(2*theta)/(4*sigx2) - np.sin(2*theta)/(4*sigy2)
	
	expo = -a*(x-x0)**2 - b*(y-y0)**2 - 2*c*(x-x0)*(y-y0)
	return A*np.exp(expo)    

def getest2DGF(x, y, I):
	"""Given the meshgrid (Xg,Yg) in x and y and the image intensities on that grid in I
	then use the moments of your image to calculate constants a, b and c from which we can
	calculate the semi major axis length (1 sigma) and semi minor axis legth (1 sigma) and
	the angle between the major axis and the positive x axis measured counter clockwise.
	If values a, b and c do not comply to the conditions for an ellipse, then return None.
	The calling environment should check the return value."""
	
	M0 = I.sum()
	x0 = (x*I).sum()/M0
	y0 = (y*I).sum()/M0
	Mxx = (x*x*I).sum()/M0 - x0*x0
	Myy = (y*y*I).sum()/M0 - y0*y0
	Mxy = (x*y*I).sum()/M0 - x0*y0
	D = 2*(Mxx*Myy-Mxy*Mxy)
	a = Myy/D
	b = Mxx/D
	c = -Mxy/D
	if a*b-c*c < 0 or a <= 0 or b <= 0:
		return None
	
	# Find the area of one pixel expressed in grids to find amplitude A
	Nx = x[0].size
	Ny = y[:,0].size
	dx = abs(x[0,0]-x[0,-1])/Nx
	dy = abs(y[0,0]-y[-1,-1])/Ny
	A = dx*dy*M0*(a*b-c*c)**0.5/np.pi
	
	p = ((a-b)**2+4*c*c)**0.5   
	theta = np.degrees(0.5*np.arctan(2*c/(a-b))) 
	if a-b > 0: # Not HW1 but the largest axis corresponds to theta.
		theta += 90.0
	if theta < 0:
		theta += 180
	
	# Major and minor axis lengths
	major = (2/(a+b-p))**0.5
	minor = (2/(a+b+p))**0.5
	sx = major
	sy = minor
	return A, x0, y0, sx, sy, theta

files = sorted(glob.glob("*.rescaled"))
a = psrchive.Archive_load(files[0])
a.remove_baseline()
a.tscrunch()
a.pscrunch()
a.centre()
data2 = a.get_data()
nsub, npol, nchan, nbin = data2.shape
# ms per bin
period = a.integration_length()
mspb = 1000*period/nbin
# folded profile and peak index
F = np.mean(data2[0,0,:,:], axis=0)
peak_idx = np.argmax(F)
if sys.argv[1] == 'pulse_65080037.calib.rescaled':
		peaks, _ = find_peaks(F)
		peak_idx = np.where(F==np.sort(F[peaks])[-2])[0][0]
# phase window
width = np.round(8*peak_widths(F, np.array([peak_idx]), rel_height=0.5)[0]).astype(int)
# on-pulse phase bin start and finish
ps = int(np.round(peak_idx - width/2))
pf = int(np.round(peak_idx + width/2))
### FREQ ZOOM
bw = a.get_bandwidth()
cf = a.get_centre_frequency()
#lowest observed frequency
min_freq = cf-bw/2
# fscrunching factor
f1 = float(sys.argv[3])
f2 = float(sys.argv[4])
f_scr = bw/a.get_nchan()
fs = int((f1-min_freq)/f_scr)
ff = int((f2-min_freq)/f_scr)
# intensity
P = data2[0,0,fs:ff,ps:pf]
nchan, nbin = P.shape

if os.path.isfile('duration_drift.npy'):
	duration_drift = np.load('duration_drift.npy')
else:
	duration_drift = []
	for f in sorted(glob.glob("corr_2D_%s.npy"%sys.argv[1].split(os.extsep, 1)[0])):
		corr_2D = np.load('corr_2D_%s.npy'%sys.argv[1].split(os.extsep, 1)[0])
		corr_2D[nchan-18:nchan+18,:] = (corr_2D[nchan-19,:] + corr_2D[nchan+19,:])/2
		sum_corr_2D_freq = np.mean(corr_2D, axis=1)
		#print(sum_corr_2D_freq.shape)
		#print(np.amax(sum_corr_2D_freq), np.argmax(sum_corr_2D_freq))
		# summed phase auto-correlation
		sum_corr_2D_time = np.mean(corr_2D, axis=0)
		## fit 2D gaussian for drift-rate
		x = np.linspace(0, corr_2D.shape[1], corr_2D.shape[1])
		y = np.linspace(0, corr_2D.shape[0], corr_2D.shape[0])
		Xg, Yg = np.meshgrid(x, y)
		# estimate initial parameters
		Ae, x0e, y0e, sxe, sye, thetae = getest2DGF(Xg, Yg, corr_2D)
		#print("Found initial 2D gaussian estimates: ", Ae, x0e, y0e, sxe, sye, thetae)
		# estimated 2D gaussian
		# fit 2D gaussian with lmfit
		fmodel = Model(RotGauss2D, independent_vars=('x','y'))
		result = fmodel.fit(corr_2D, x=Xg, y=Yg, A=Ae, x0=x0e, y0=y0e, sigma_x=sxe, sigma_y=sye, theta=thetae)
		#print(lmfit.report_fit(result))
		corr_2D_model = RotGauss2D(Xg, Yg, result.best_values['A'], result.best_values['x0'], result.best_values['y0'], 
					   result.best_values['sigma_x'], result.best_values['sigma_y'], result.best_values['theta'])
		FWHM = 2*np.sqrt(np.log(2))
		#ell = Ellipse(xy=(result.best_values['x0'], result.best_values['y0']), width=FWHM*result.best_values['sigma_x'], 
		#	      height=FWHM*result.best_values['sigma_y'], angle=result.best_values['theta'], fill=False, color='r')
		# drift rate from fit
		theta = result.best_values['theta']
		if 180 > theta > 90:
			drift_rate = np.round(-np.tan((180-theta)*np.pi/180)*(bw/nchan)/(mspb), 2)
		if 0 < theta < 90:
			drift_rate = np.round(np.tan(theta*np.pi/180)*(bw/nchan)/(mspb), 2)
		if theta == 90:
			drift_rate = 'No Drift'
		if theta == 0:
			drift_rate = 0

		sigma_bw_tot = np.abs(gauss_fit(x,sum_corr_2D_freq)[-1])
		#append duration and drift rate
		duration_drift.append(np.array([np.round(FWHM*sigma_bw_tot*(bw/nchan), 1), drift_rate]))
	np.save('duration_drift.npy', duration_drift)

### PLOTTING ###
A4x, A4y = 8.27, 11.69
fig = plt.figure(figsize=(A4y, A4x), dpi=600)

for i in range(len(files)):
	plt.plot(duration_drift[i][0], duration_drift[i][1], 'o', color='k', markersize=2)
	plt.xlabel('Duration (ms)')
	plt.ylabel('Drift Rate (MHz/ms)')
plt.savefig('duration_drift.png', dpi=600)
print('duration_drift.png')