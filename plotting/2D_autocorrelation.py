import os
import sys

import lmfit
from lmfit.model import save_modelresult
import matplotlib.pyplot as plt
import numpy as np
import psrchive
from uncertainties import unumpy as unp
from uncertainties import ufloat
from lmfit.lineshapes import lorentzian
from lmfit import Model
from matplotlib.patches import Ellipse
from matplotlib import gridspec
from scipy import signal
from scipy.optimize import curve_fit
from scipy.signal import peak_widths, find_peaks
from mpl_toolkits.mplot3d import Axes3D

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

def RotGauss2Dnew(x, y, A, x0, y0, sigma_x, sigma_y, theta):
	"""Rotated 2D Gaussian function"""
	theta = np.radians(theta)
	k = sigma_x/sigma_y - sigma_y/sigma_x
	p = (k*np.sin(2*theta)/2)/np.sqrt(1+(k*np.sin(2*theta)/2)**2)
	w_t = sigma_x*sigma_y*np.sqrt((1+(k*np.sin(2*theta)/2)**2)/(sigma_x**2*np.sin(theta)**2+sigma_y**2*np.cos(theta)**2))
	sigma_v = sigma_x*sigma_y*np.sqrt((1+(k*np.sin(2*theta)/2)**2)/(sigma_x**2*np.cos(theta)**2+sigma_y**2*np.sin(theta)**2))
	
	expo = -1/(1-p**2)*(((x-x0)**2)/(2*w_t**2) - p*(x-x0)*(y-y0)/(w_t*sigma_v) + ((y-y0)**2)/(2*sigma_v**2))
	return A*np.exp(expo)    

def Gauss2D_dt(x, y, A, x0, y0, sigma_x, sigma_y, d_t):
	"""2D Gaussian function with time delay for intra-burst drift"""
	a = 1/(2*sigma_x**2)
	b = -d_t/(2*sigma_x**2)
	c = d_t**2/(2*sigma_x**2) + 1/(2*sigma_y**2)

	#expo = -(a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2)
	expo = -(x-x0-d_t*(y-y0))**2/(2*sigma_x**2) - (y-y0)**2/(2*sigma_y**2)
	return A*np.exp(expo)

def Gauss2D_dv(x, y, A, x0, y0, sigma_x, sigma_y, d_t):
	"""2D Gaussian function with time delay for sub-burst drift"""
	d_v = d_t/(d_t**2+sigma_x**2/sigma_y**2)
	w_v = sigma_x/np.sqrt(1+(d_t*sigma_y/sigma_x)**2)
	w_t = np.sqrt(sigma_x**2+(d_t*sigma_y)**2)

	a = 1/(2*w_t**2) + d_v**2/(2*w_v**2)
	b = -d_v/(2*w_v**2)
	c = 1/(2*w_v**2)

	expo = -(a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2)
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


f1 = float(sys.argv[3])
f2 = float(sys.argv[4])
p = sys.argv[2]
if p == "I":
	a = psrchive.Archive_load(sys.argv[1])
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
	if sys.argv[1] == 'pulse_65080037.calib.rescaled':
			peaks, _ = find_peaks(F)
			peak_idx = np.where(F==np.sort(F[peaks])[-2])[0][0]

	# phase window
	width = np.round(8*peak_widths(F, np.array([peak_idx]), rel_height=0.5)[0]).astype(int)

	# on-pulse phase bin start and finish
	ps = int(np.round(peak_idx - width/2))
	pf = int(np.round(peak_idx + width/2))

	h=0.3*np.amax(F[ps:pf])
	peaks, _ = find_peaks(F[ps:pf], height=h)
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
	### FREQ ZOOM
	bw = a.get_bandwidth()
	cf = a.get_centre_frequency()
	#lowest observed frequency
	min_freq = cf-bw/2
	# fscrunching factor
	f_scr = bw/a.get_nchan()

	fs = int((f1-min_freq)/f_scr)
	ff = int((f2-min_freq)/f_scr)

	# intensity
	P = data2[0,0,fs:ff,ps:pf]

else:
	a = psrchive.Archive_load(sys.argv[1])
	a.remove_baseline()
	a.tscrunch()
	a.centre()
	data2 = a.get_data()
	nsub, npol, nchan, nbin = data2.shape

	# folded profile and peak index
	F = data2.mean(axis=(1,2))[0]
	peak_idx = np.argmax(F)

	# on-pulse phase start and finish
	width = np.round(6*peak_widths(F[peak_idx], peak_idx, rel_height=0.5)[0]).astype(int)

	# phase window
	ps = int(np.round(peak_idx - width/2))
	pf = int(np.round(peak_idx + width/2))

	### FREQ ZOOM
	bw = a.get_bandwidth()
	cf = a.get_centre_frequency()
	f_scr = bw/a.get_nchan()

	fs = int((f1-(cf-bw/2))/f_scr)
	ff = int((f2-(cf-bw/2))/f_scr)

	# polarisations
	if p == "SI":
		P = data2[0,0,fs:ff,ps:pf]
	if p == "SQ":
		P = data2[0,1,fs:ff,ps:pf]
	if p == "SU":
		P = data2[0,2,fs:ff,ps:pf]
	if p == "L":
		P = np.sqrt(data2[0,1,fs:ff,ps:pf]**2+data2[0,2,fs:ff,ps:pf]**2)
	if p == "SV":
		P = data2[0,3,fs:ff,ps:pf]
nchan, nbin = P.shape

## 1D auto-correlation of each frequency channel
if os.path.isfile('freq_corr_%s.npy'%sys.argv[1].split(os.extsep, 1)[0]):
	freq_corr = np.load('freq_corr_%s.npy'%sys.argv[1].split(os.extsep, 1)[0])
else:
	freq_corr = []
	for i in range(nchan):
		freq_corr.append(signal.correlate(P[i,:], P[i,:], mode='full', method='direct'))
	np.save("freq_corr_%s.npy"%sys.argv[1].split(os.extsep, 1)[0], freq_corr)
# summed frequency auto-correlation
sum_freq_corr_freq = np.mean(freq_corr, axis=1)
sum_freq_corr_time = np.mean(freq_corr, axis=0)
# fit 1D gaussian to summed-frequency auto-correlation

## 1D auto-correlation of each phase bin
if os.path.isfile('phase_corr_%s.npy'%sys.argv[1].split(os.extsep, 1)[0]):
	phase_corr = np.load('phase_corr_%s.npy'%sys.argv[1].split(os.extsep, 1)[0])
else:
	phase_corr = []
	for i in range(nbin):
		phase_corr.append(signal.correlate(P[:,i], P[:,i], mode='full', method='direct'))
	np.save("phase_corr_%s.npy"%sys.argv[1].split(os.extsep, 1)[0], phase_corr)
phase_corr = np.rot90(phase_corr)
phase_corr[3310:3344,:] = (phase_corr[3309,:] + phase_corr[3345,:])/2
# summed phase auto-correlation
sum_phase_corr_freq = np.mean(phase_corr, axis=1)
sum_phase_corr_time = np.mean(phase_corr, axis=0)


## 2D auto-correlation
if os.path.isfile('corr_2D_%s.npy'%sys.argv[1].split(os.extsep, 1)[0]):
	corr_2D = np.load('corr_2D_%s.npy'%sys.argv[1].split(os.extsep, 1)[0])
else:	
	corr_2D = signal.correlate2d(P, P, mode='full', boundary='fill', fillvalue=0)
	np.save("corr_2D_%s.npy"%sys.argv[1].split(os.extsep, 1)[0], corr_2D)
corr_2D[nchan-18:nchan+18,:] = (corr_2D[nchan-19,:] + corr_2D[nchan+19,:])/2
sum_corr_2D_freq = np.mean(corr_2D, axis=1)
#print(sum_corr_2D_freq.shape)
#print(np.amax(sum_corr_2D_freq), np.argmax(sum_corr_2D_freq))
# summed phase auto-correlation
sum_corr_2D_time = np.mean(corr_2D, axis=0)
## fit 2D gaussian for drift-rate
if len(peaks) <= 1:
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
	with open('lmfit_result.txt', 'w') as fh:
		fh.write(result.fit_report())
	# pull theta error from file
	line = open('lmfit_result.txt').readlines()[18]
	theta_error = float(line[line.index("-"):line.index("(")][1:])*np.pi/180

	theta = result.best_values['theta']*np.pi/180
	if 180 > theta > 0:
		dr = 1/unp.tan(np.array([ufloat(theta,theta_error)]))*mspb/(bw/nchan)
		drift_rate = float(str(dr)[1:-1].split("+/-")[0])
	if theta == 90:
		print("No Drift")
		sys.exit()
	if theta == 0:
		print("drift rate = 0")
	#drift_err = np.abs(theta_error*(np.pi**2/180**2)*(1/np.cos(theta*np.pi/180))**2)*(bw/nchan)/(mspb)
	drift_err = float(str(dr)[1:-1].split("+/-")[1])
else:
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
	#ell = Ellipse(xy=(result.best_values['x0'], result.best_values['y0']), width=FWHM*result.best_values['sigma_x'], 
	#	      height=FWHM*result.best_values['sigma_y'], angle=result.best_values['theta'], fill=False, color='r')
	# drift rate from fit
	# save fit report to a file:
	with open('lmfit_result.txt', 'w') as fh:
		fh.write(result.fit_report())
	# pull theta error from file
	line = open('lmfit_result.txt').readlines()[18]
	theta_error = float(line[line.index("-"):line.index("(")][1:])*np.pi/180

	theta = result.best_values['theta']*np.pi/180
	if 180 > theta > 0:
		dr = unp.tan(np.array([ufloat(theta,theta_error)]))*(bw/nchan)/mspb
		drift_rate = np.round(float(str(dr)[1:-1].split("+/-")[0]), 2)
	if theta == 90:
		print("No Drift")
		sys.exit()
	if theta == 0:
		print("drift rate = 0")
		sys.exit()
	#drift_err = np.abs(-theta_error*(np.pi**2/180**2)*(1/np.sin(theta*np.pi/180))**2)*(bw/nchan)/(mspb)
	drift_err = np.round(float(str(dr)[1:-1].split("+/-")[1]), 2)



FWHM = 2*np.sqrt(np.log(2))
#if len(peaks) == 1:
#	if 180 > theta > 90:
#		drift_rate = np.round(-1/((np.tan((180-theta))*np.pi/180)*(bw/nchan)/(mspb)), 2)
#	if 0 < theta < 90:
#		drift_rate = np.round(1/((np.tan(theta*np.pi/180))*(bw/nchan)/(mspb)), 2)
#	if theta == 90:
#		print("No Drift")
#		sys.exit()
#	if theta == 0:
#		print("drift rate = 0")
#		sys.exit()
#	### MODEL FROM JAHNS ET AL. 2023 ###
#	fmodel = Model(Gauss2D_dt, independent_vars=('x','y'))
#	result = fmodel.fit(corr_2D, x=Xg*mspb, y=Yg+min_freq, A=Ae, x0=x0e*mspb, y0=y0e+min_freq, sigma_x=sye+min_freq, sigma_y=sxe*mspb, d_t=drift_rate)
#	print(lmfit.report_fit(result))
#	corr_2D_model = Gauss2D_dt(Xg*mspb, Yg+min_freq, result.best_values['A'], result.best_values['x0'], result.best_values['y0'], 
#				   result.best_values['sigma_x'], result.best_values['sigma_y'], result.best_values['d_t'])
#	drift_rate = np.round(result.best_values['d_t'], 2)
#else:
#	if 180 > theta > 90:
#		drift_rate = np.round(-1/((np.tan((180-theta))*np.pi/180)*(bw/nchan)/(mspb)), 2)
#	if 0 < theta < 90:
#		drift_rate = np.round(1/((np.tan(theta*np.pi/180))*(bw/nchan)/(mspb)), 2)
#	if theta == 90:
#		drift_rate = 'No Drift'
#	if theta == 0:
#		drift_rate = 0
#	### MODEL FROM JAHNS ET AL. 2023 ###
#	fmodel = Model(Gauss2D_dv, independent_vars=('x','y'))
#	result = fmodel.fit(corr_2D, x=Xg, y=Yg, A=Ae, x0=x0e*mspb, y0=y0e, sigma_x=sxe*mspb, sigma_y=sye, d_t=drift_rate)
#	#print(lmfit.report_fit(result))
#	corr_2D_model = Gauss2D_dv(Xg, Yg, result.best_values['A'], result.best_values['x0'], result.best_values['y0'], 
#				   result.best_values['sigma_x'], result.best_values['sigma_y'], result.best_values['d_t'])
#	drift_rate = np.round(result.best_values['d_t'], 2)






# 3D plot of 2D auto-correlation
#fig = plt.figure()
#frame = fig.add_subplot(1,1,1, projection='3d', azim=-45, elev=30)
#step = 21
#frame.plot_surface(Xg, Yg, corr_2D_model, cmap='jet')
#frame.set_xlabel('X axis')
#frame.set_ylabel('Y axis')
#frame.set_zlabel('Z axis')
#plt.show()


### PLOTTING ###
A4x, A4y = 8.27, 11.69
fig = plt.figure(figsize=(A4x, A4x), dpi=600)
g = gridspec.GridSpec(ncols=3, nrows=4, hspace=0., wspace=0., 
			  height_ratios=[0.25,0.25,1,1], width_ratios=[1,1,0.5])

# plot dynamic spectrum
ms_tick = nbin*mspb
dy_spec_xticks = np.round(np.linspace(0,ms_tick,num=5),2)
dy_spec_xticks_x = np.linspace(0,pf-ps-1,num=len(dy_spec_xticks))
dy_spec_yticks = np.linspace(f1,f2, num=7).astype(int)
dy_spec_yticks_y = np.linspace(0,ff-fs-1, len(dy_spec_yticks))

ax_3_0 = fig.add_subplot(g[3,0])
ax_3_0.imshow(P, cmap='Greys', aspect='auto', origin='lower', interpolation='none', vmax=0.9*np.amax(P))
ax_3_0.set_xticks(dy_spec_xticks_x)
ax_3_0.set_xticklabels(dy_spec_xticks)
ax_3_0.set_yticks(dy_spec_yticks_y)
ax_3_0.set_yticklabels(dy_spec_yticks)
ax_3_0.set_xlabel('Time (ms)')
ax_3_0.set_ylabel('Frequency (MHz)')

# plot frequency autocorrelation
f_cor_xticks = np.round(np.linspace(-ms_tick, ms_tick, num=5),2)
f_cor_xticks_x = np.linspace(0,2*(nbin-1),num=len(f_cor_xticks))

ax_3_1 = fig.add_subplot(g[3,1])
ax_3_1.imshow(freq_corr, cmap='Reds', aspect='auto', origin='lower', interpolation='none')
ax_3_1.set_xticks(f_cor_xticks_x[1:-1])
ax_3_1.set_xticklabels(f_cor_xticks[1:-1])
ax_3_1.set_yticklabels([])
ax_3_1.set_yticks(dy_spec_yticks_y)	
ax_3_1.set_xlabel('Time Lag (ms)')

# plot phase autocorrelation
p_cor_yticks = np.linspace(-(f2-f1-1), f2-f1-1, num=7).astype(int)
p_cor_yticks_y = np.linspace(0,2*(nchan-1), len(p_cor_yticks))

ax_2_0 = fig.add_subplot(g[2,0])
ax_2_0.imshow(phase_corr, cmap='Blues', aspect='auto', origin='lower', interpolation='none')
ax_2_0.set_xticklabels([])
ax_2_0.set_xticks(dy_spec_xticks_x)
ax_2_0.set_yticks(p_cor_yticks_y[1:-1])
ax_2_0.set_yticklabels(p_cor_yticks[1:-1])
ax_2_0.set_ylabel('Frequency Lag (MHz)')

# plot 2D auto-correlation
ax_2_1 = fig.add_subplot(g[2,1])
ax_2_1.imshow(corr_2D, cmap='Purples', aspect='auto', origin='lower', interpolation='none')
# contour plot
peak = np.amax(corr_2D_model)
ax_2_1.contour(Xg,Yg, corr_2D_model, levels=[peak/2], colors='k', linewidths=1)
# drift-rate text box
props = dict(boxstyle='square', facecolor='white', alpha=0.5)
if len(peaks) <= 1:
	ax_2_1.text(0.5, 0.95, r'$\frac{d\nu}{dt}$ = %s $\times 10^{-5}$ ms MHz$^{-1}$' % np.round(drift_rate*10**5, 2), 
		transform=ax_2_1.transAxes, fontsize=10, verticalalignment='top', horizontalalignment ='center', 
		bbox=props, family='serif', fontweight='ultralight')
else:
	ax_2_1.text(0.5, 0.95, r'$\frac{d\nu}{dt}$ = %s MHz ms$^{-1}$' % drift_rate, 
		transform=ax_2_1.transAxes, fontsize=10, verticalalignment='top', horizontalalignment ='center', 
		bbox=props, family='serif', fontweight='ultralight')
#ax_2_1.add_patch(ell)
ax_2_1.set_xticklabels([])
ax_2_1.set_xticks(f_cor_xticks_x[1:-1])
ax_2_1.set_yticklabels([])
ax_2_1.set_yticks(p_cor_yticks_y[1:-1])


## Summed auto-correlations
# plot summed frequency autocorrelation in frequency
ax_3_2 = fig.add_subplot(g[3,2])
ax_3_2.step(sum_freq_corr_freq, np.arange(len(sum_freq_corr_freq)), 
	color='red', where='mid', lw=0.5)
ax_3_2.margins(y=0)
ax_3_2.set_xticklabels([])
ax_3_2.set_yticklabels([])
ax_3_2.set_yticks(dy_spec_yticks_y)
ax_3_2.set_xlabel('(arb. units)')

# plot summed frequency autocorrelation in time
ax_1_1 = fig.add_subplot(g[1,1])
ax_1_1.step(np.arange(len(sum_freq_corr_time)), sum_freq_corr_time, 
	color='red', where='mid', lw=0.5)
# plot 1D gaussian fit
x = np.arange(len(sum_freq_corr_time))
ax_1_1.plot(x, gauss(x, *gauss_fit(x, sum_freq_corr_time)), color='red', lw=1)
sigma_dur_sub = np.abs(gauss_fit(x,sum_freq_corr_time)[-1])
ax_1_1.text(0.05, 0.95, r'$\Delta t_{{\mathrm{sub}}}$'
	    				'\n'
						'%s ms'%np.round(FWHM*sigma_dur_sub*mspb, 1), 
						transform=ax_1_1.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left', 
						color='red', family='serif', fontweight='ultralight')
ax_1_1.margins(x=0)
ax_1_1.set_xticklabels([])
ax_1_1.set_xticks(f_cor_xticks_x[1:-1])
ax_1_1.set_yticklabels([])

# plot summed phase autocorrelation in frequency
ax_2_2 = fig.add_subplot(g[2,2])
ax_2_2.step(sum_phase_corr_freq, np.arange(len(sum_phase_corr_freq)), 
	color='blue', where='mid', lw=0.5)
# plot 1D gaussian fit
x = np.arange(len(sum_phase_corr_freq))
ax_2_2.plot(gauss(x, *gauss_fit(x, sum_phase_corr_freq)) , x, color='blue', lw=1)
sigma_bw_sub = np.abs(gauss_fit(x,sum_phase_corr_freq)[-1])
ax_2_2.text(0.95, 0.05, r'$\Delta\nu_{{\mathrm{sub}}}$'
	    				'\n'
						 '%s MHz'%np.round(FWHM*sigma_bw_sub*(bw/nchan), 1), 
						 transform=ax_2_2.transAxes, fontsize=10, verticalalignment='bottom', horizontalalignment='right', 
						 color='blue', family='serif', fontweight='ultralight')
ax_2_2.margins(y=0)
ax_2_2.set_xticklabels([])
ax_2_2.set_yticklabels([])
ax_2_2.set_yticks(p_cor_yticks_y[1:-1])

# plot summed phase autocorrelation in time
ax_02_0 = fig.add_subplot(g[0:2,0])
ax_02_0.step(np.arange(len(sum_phase_corr_time)), sum_phase_corr_time, 
	color='blue', where='mid', lw=0.5)
ax_02_0.margins(x=0)
ax_02_0.set_xticklabels([])
ax_02_0.set_xticks(dy_spec_xticks_x)
ax_02_0.set_yticklabels([])
ax_02_0.set_ylabel('(arb. units)')

# plot summed 2D autocorrelation in frequency
ax_2_2.step(sum_corr_2D_freq, np.arange(len(sum_corr_2D_freq)), 
	color='purple', where='mid', lw=0.5)
# plot 1D gaussian fit
x = np.arange(len(sum_corr_2D_freq))
ax_2_2.plot(gauss(x, *gauss_fit(x, sum_corr_2D_freq)), x, color='purple', lw=1)
sigma_bw_tot = np.abs(gauss_fit(x,sum_corr_2D_freq)[-1])
ax_2_2.text(0.95, 0.95, r'$\Delta\nu_{{\mathrm{tot}}}$'
	    				'\n'
						 '%s MHz'%np.round(FWHM*sigma_bw_tot*(bw/nchan), 1), 
						 transform=ax_2_2.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right', 
						 color='purple', family='serif', fontweight='ultralight')

# plot summed 2D autocorrelation in time
ax_0_1 = fig.add_subplot(g[0,1])
ax_0_1.step(np.arange(len(sum_corr_2D_time)), sum_corr_2D_time, 
	color='purple', where='mid', lw=0.5)
ax_0_1.margins(x=0)
ax_0_1.set_xticklabels([])
ax_0_1.set_xticks(f_cor_xticks_x[1:-1])
ax_0_1.set_yticklabels([])
# plot 1D gaussian fit
x = np.arange(len(sum_corr_2D_time))
ax_0_1.plot(x, gauss(x, *gauss_fit(x, sum_corr_2D_time)), color='purple', lw=1)
sigma_dur_tot = np.abs(gauss_fit(x,sum_corr_2D_time)[-1])
ax_0_1.text(0.05, 0.95, r'$\Delta t_{{\mathrm{tot}}}$'
	    				'\n'
						 '%s ms'%np.round(FWHM*sigma_dur_tot*mspb, 1), 
						 transform=ax_0_1.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left', 
						 color='purple', family='serif', fontweight='ultralight')

# PRINT DRIFT RATE BW AND DURATION
if len(peaks) <= 1:
	print("DRIFT = %s \pm %s ms/MHz"%(drift_rate, drift_err))
else:
	print("DRIFT = %s \pm %s MHz/ms"%(drift_rate, drift_err))
print("BW_tot = %s \pm %s MHz" %(np.round(FWHM*sigma_bw_tot*(bw/nchan), 1), np.round((sigma_bw_tot/np.sqrt(len(sum_corr_2D_freq)))*(bw/nchan), 2)))
print("BW_sub = %s \pm %s MHz" %(np.round(FWHM*sigma_bw_sub*(bw/nchan), 1), np.round((sigma_bw_sub/np.sqrt(len(sum_phase_corr_freq)))*(bw/nchan), 2)))
print("Dur_tot = %s \pm %s ms" %(np.round(FWHM*sigma_dur_tot*mspb, 1), np.round((sigma_dur_tot/np.sqrt(len(sum_corr_2D_time)))*mspb, 2)))
print("Dur_sub = %s \pm %s ms" %(np.round(FWHM*sigma_dur_sub*mspb, 1), np.round((sigma_dur_sub/np.sqrt(len(sum_freq_corr_time)))*mspb, 2)))
plt.savefig(sys.argv[0].split(os.extsep, 1)[0]+'_%s_'%sys.argv[2]+sys.argv[1].split(os.extsep, 1)[0]+'.png', dpi=600, bbox_inches='tight')
print(sys.argv[0].split(os.extsep, 1)[0]+'_%s_'%sys.argv[2]+sys.argv[1].split(os.extsep, 1)[0]+'.pdf')

print("$%s \\times 10^{-5} \pm %s \\times 10^{-5} $ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ \\\\"%(np.round(drift_rate*1e5, 2), np.round(drift_err*1e5, 3), 
												np.round(FWHM*sigma_bw_tot*(bw/nchan), 1), np.round((sigma_bw_tot/np.sqrt(len(sum_corr_2D_freq)))*(bw/nchan), 2),
												np.round(FWHM*sigma_bw_sub*(bw/nchan), 1), np.round((sigma_bw_sub/np.sqrt(len(sum_phase_corr_freq)))*(bw/nchan), 2),
												np.round(FWHM*sigma_dur_tot*mspb, 1), np.round((sigma_dur_tot/np.sqrt(len(sum_corr_2D_time)))*mspb, 2),
												np.round(FWHM*sigma_dur_sub*mspb, 1), np.round((sigma_dur_sub/np.sqrt(len(sum_freq_corr_time)))*mspb, 2))
												)
