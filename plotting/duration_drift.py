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



if os.path.isfile('duration_intra_drift.npy') and os.path.isfile('duration_sub_drift.npy'):
	duration_intra_drift = np.load('duration_intra_drift.npy')
	duration_sub_drift = np.load('duration_sub_drift.npy')
	phase = np.load('phase.npy')
	duration_intra_drift = np.array(duration_intra_drift)
	duration_sub_drift = np.array(duration_sub_drift)
	phase = np.array(phase)
else:
	duration_intra_drift = []
	duration_sub_drift = []
	phase = []
	files = sorted(glob.glob("*.rescaled"))
	for f in files:
		a = psrchive.Archive_load(f)
		a.remove_baseline()
		a.tscrunch()
		a.pscrunch()
		data2 = a.get_data()
		nsub, npol, nchan, nbin = data2.shape
		# ms per bin
		period = a.integration_length()
		mspb = 1000*period/nbin
		# folded profile and peak index
		F = np.mean(data2[0,0,:,:], axis=0)
		peak_idx = np.argmax(F)
		# phase window
		width = np.round(8*peak_widths(F, np.array([peak_idx]), rel_height=0.5)[0]).astype(int)
		# on-pulse phase bin start and finish
		ps = int(np.round(peak_idx - width/2))
		pf = int(np.round(peak_idx + width/2))
		phase.append(np.array([ps, pf]))

		bw = a.get_bandwidth()
		# intensity
		P = data2[0,0,:,ps:pf]
		nchan, nbin = P.shape

		corr_2D = np.load('corr_2D_'+f.split('.')[0]+'.npy')
		corr_2D[nchan-18:nchan+18,:] = (corr_2D[nchan-19,:] + corr_2D[nchan+19,:])/2
		#sum_corr_2D_freq = np.mean(corr_2D, axis=1)

		# summed phase auto-correlation
		sum_corr_2D_time = np.mean(corr_2D, axis=0)
		## fit 2D gaussian for drift-rate
		x = np.linspace(0, corr_2D.shape[1], corr_2D.shape[1])
		y = np.linspace(0, corr_2D.shape[0], corr_2D.shape[0])
		Xg, Yg = np.meshgrid(x, y)

		# estimate initial parameters
		if getest2DGF(Xg, Yg, corr_2D) != None:
			Ae, x0e, y0e, sxe, sye, thetae = getest2DGF(Xg, Yg, corr_2D)

			# estimated 2D gaussian
			# fit 2D gaussian with lmfit
			fmodel = Model(RotGauss2D, independent_vars=('x','y'))
			result = fmodel.fit(corr_2D, x=Xg, y=Yg, A=Ae, x0=x0e, y0=y0e, sigma_x=sxe, sigma_y=sye, theta=thetae)

			corr_2D_model = RotGauss2D(Xg, Yg, result.best_values['A'], result.best_values['x0'], result.best_values['y0'], 
						   result.best_values['sigma_x'], result.best_values['sigma_y'], result.best_values['theta'])
			FWHM = 2*np.sqrt(np.log(2))

			theta = result.best_values['theta']
			h=0.3*np.amax(F[ps:pf])
			peaks, _ = find_peaks(F[ps:pf], height=h)
			if len(peaks) == 1:
				if 180 > theta > 90:
					drift_rate = np.round(1/((np.tan((180-theta))*np.pi/180)*(bw/nchan)/(mspb)), 2)
				if 0 < theta < 90:
					drift_rate = np.round(1/((np.tan(theta*np.pi/180))*(bw/nchan)/(mspb)), 2)
				if theta == 90:
					print("No Drift")
				if theta == 0:
					print("drift rate = 0")

			else:
				if 180 > theta > 90:
					drift_rate = np.round(-np.tan((180-theta)*np.pi/180)*(bw/nchan)/(mspb), 2)
				if 0 < theta < 90:
					drift_rate = np.round(np.tan(theta*np.pi/180)*(bw/nchan)/(mspb), 2)
				if theta == 90:
					print("No Drift")
				if theta == 0:
					print("drift rate = 0")

			x = np.arange(len(sum_corr_2D_time))
			sigma_dur_tot = np.abs(gauss_fit(x,sum_corr_2D_time)[-1])
			Dur_tot = np.round(FWHM*sigma_dur_tot*mspb, 3)
			if -26000 < drift_rate < 26000:
				# peak minimas
				mins, _ = find_peaks(-F[ps:pf])
				# associate peaks with minimas
				try:
					if peaks[-1] > mins[-1]:
						peaks = peaks[:-1]
					if peaks[0] < mins[0]:
						peaks = peaks[1:]
				except IndexError:
					pass
				if (len(peaks) == 1):
					#append duration and drift rate
					duration_intra_drift.append(np.array([Dur_tot, drift_rate, int(f.split('_')[-1].split('.')[0])]))
				if (len(peaks) > 1 and Dur_tot < 10):
					#append duration and drift rate
					duration_sub_drift.append(np.array([Dur_tot, drift_rate, int(f.split('_')[-1].split('.')[0])]))
				#append duration and drift rate
	np.save('duration_intra_drift.npy', duration_intra_drift)
	np.save('duration_sub_drift.npy', duration_sub_drift)
	np.save('phase.npy', phase)
	duration_intra_drift = np.array(duration_intra_drift)
	duration_sub_drift = np.array(duration_sub_drift)
	phase = np.array(phase)


### PLOTTING ###
A4x, A4y = 8.27, 11.69
if len(duration_intra_drift) != 0:
	fig = plt.figure(figsize=(A4y, A4x), dpi=600)
	plt.plot(duration_intra_drift[:,0], duration_intra_drift[:,1], 'o', color='k', markersize=2)
	#plt.xscale('log')
	#plt.yscale('log')
	#plt.ylim(-10, 10)
	plt.xlabel('Duration (ms)')
	plt.ylabel('Drift Rate (ms/MHz)')
	plt.savefig('duration_intra_drift.png', dpi=600, bbox_inches='tight')
	print('duration_intra_drift.png')

if len(duration_sub_drift) != 0:
	fig = plt.figure(figsize=(A4y, A4x), dpi=600)
	plt.plot(duration_sub_drift[:,0], duration_sub_drift[:,1], 'o', color='k', markersize=2)
	#plt.xscale('log')
	#plt.yscale('log')
	plt.xlabel('Duration (ms)')
	plt.ylabel('Drift Rate (MHz/ms)')
	plt.savefig('duration_sub_drift.png', dpi=600, bbox_inches='tight')
	print('duration_sub_drift.png')