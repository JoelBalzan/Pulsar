import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import glob
import os
from matplotlib import gridspec


profiles_P970 = np.load('pulse_profiles_P970.npy')
nfile_P970, nbin_P970 = np.shape(profiles_P970)
on_s_P970 = int(0.53*nbin_P970)
on_f_P970 = int(0.78*nbin_P970)
mean_flux_P970 = []
for i in range(nfile_P970):
	flux = np.mean(profiles_P970[i][on_s_P970:on_f_P970], axis=0)
	mean_flux_P970.append(flux)
avg_P970 = np.mean(mean_flux_P970)
print("Mean P970 = ", avg_P970, " Jy")
print("Peak E = ", np.max(mean_flux_P970)/avg_P970, "<E>")


profiles_PX500_38329 = np.load('pulse_profiles_PX500_38329.npy')
nfile_PX500_38329, nbin_PX500_38329 = np.shape(profiles_PX500_38329)
on_s_PX500_38329 = int(0.27*nbin_PX500_38329)
on_f_PX500_38329 = int(0.51*nbin_PX500_38329)
mean_flux_PX500_38329 = []
for i in range(nfile_PX500_38329):
	flux = np.mean(profiles_PX500_38329[i][on_s_PX500_38329:on_f_PX500_38329], axis=0)
	mean_flux_PX500_38329.append(flux)
avg_PX500_38329 = np.mean(mean_flux_PX500_38329)
print("Mean PX500_38329 = ", avg_PX500_38329, " Jy")
print("Peak E = ", np.max(mean_flux_PX500_38329)/avg_PX500_38329, "<E>")


profiles_PX500_38907 = np.load('pulse_profiles_PX500_38907.npy')
nfile_PX500_38907, nbin_PX500_38907 = np.shape(profiles_PX500_38907)
on_s_PX500_38907 = int(0.68*nbin_PX500_38907)
on_f_PX500_38907 = int(0.94*nbin_PX500_38907)
mean_flux_PX500_38907 = []
for i in range(nfile_PX500_38907):
	flux = np.mean(profiles_PX500_38907[i][on_s_PX500_38907:on_f_PX500_38907], axis=0)
	mean_flux_PX500_38907.append(flux)
avg_PX500_38907 = np.mean(mean_flux_PX500_38907)
print("Mean PX500_38907 = ", avg_PX500_38907, " Jy")
print("Peak E = ", np.max(mean_flux_PX500_38907)/avg_PX500_38907, "<E>")


profiles_PX500_39167 = np.load('pulse_profiles_PX500_39167.npy')
nfile_PX500_39167, nbin_PX500_39167 = np.shape(profiles_PX500_39167)
on_s_PX500_39167 = int(0.11*nbin_PX500_39167)
on_f_PX500_39167 = int(0.33*nbin_PX500_39167)
mean_flux_PX500_39167 = []
for i in range(nfile_PX500_39167):
	flux = np.mean(profiles_PX500_39167[i][on_s_PX500_39167:on_f_PX500_39167], axis=0)
	mean_flux_PX500_39167.append(flux)
avg_PX500_39167 = np.mean(mean_flux_PX500_39167)
print("Mean PX500_39167 = ", avg_PX500_39167, " Jy")
print("Peak E = ", np.max(mean_flux_PX500_39167)/avg_PX500_39167, "<E>")



## PLOT HISTOGRAMS
A4x, A4y = 8.27, 11.69
fontsize=10
fig = plt.figure(figsize=(A4x,A4x),dpi=300)
g = gridspec.GridSpec(ncols=2, nrows=2, hspace=0.12, wspace=0.12)
#plt.rcParams["font.family"] = "serif"
#plt.rcParams["font.weight"] = "light"
ax_0_0 = fig.add_subplot(g[0, 0])
ax_0_0.hist(mean_flux_P970, bins=int(np.sqrt(len(mean_flux_P970))), edgecolor='black', color='white')
ax_0_0.axvline(avg_P970, color='r', linewidth=1, label = '%s Jy'%np.round(avg_P970, 3))
#ax.set_xscale("log")
#ax.set_yscale("log") 
ax_0_0.margins(x=0)
#ax_0_0.set_xlabel(r'$\langle{{E}}\rangle$ (Jy)', fontsize=fontsize)
#plt.ylabel('log$_{10}$(Count)')
ax_0_0.set_ylabel('Count', fontsize=fontsize)
ax_0_0.legend(loc='center right', fontsize=fontsize, handlelength=1)
ax_0_0.text(0.95, 0.95, 'MJD 58463', transform=ax_0_0.transAxes, fontsize=fontsize, verticalalignment='top', horizontalalignment='right')
#plt.title('P970')

ax_0_1 = fig.add_subplot(g[0, 1])
ax_0_1.hist(mean_flux_PX500_38329, bins=int(np.sqrt(len(mean_flux_PX500_38329))), edgecolor='black', color='white')
ax_0_1.axvline(avg_PX500_38329, color='r', linewidth=1, label = '%s Jy'%np.round(avg_PX500_38329, 3))
ax_0_1.margins(x=0)
ax_0_1.legend(loc='center right', fontsize=fontsize, handlelength=1)
ax_0_1.text(0.95, 0.95, 'MJD 58469', transform=ax_0_1.transAxes, fontsize=fontsize, verticalalignment='top', horizontalalignment='right')

ax_1_0 = fig.add_subplot(g[1, 0])
ax_1_0.hist(mean_flux_PX500_38907, bins=int(np.sqrt(len(mean_flux_PX500_38907))), edgecolor='black', color='white')
ax_1_0.axvline(avg_PX500_38907, color='r', linewidth=1, label = '%s Jy'%np.round(avg_PX500_38907, 3))
ax_1_0.margins(x=0)
ax_1_0.set_xlabel(r'$\langle{{E}}\rangle$ (Jy)', fontsize=fontsize)
ax_1_0.set_ylabel('Count', fontsize=fontsize)
ax_1_0.legend(loc='center right', fontsize=fontsize, handlelength=1)
ax_1_0.text(0.95, 0.95, 'MJD 58530', transform=ax_1_0.transAxes, fontsize=fontsize, verticalalignment='top', horizontalalignment='right')

ax_1_1 = fig.add_subplot(g[1, 1])
ax_1_1.hist(mean_flux_PX500_39167, bins=int(np.sqrt(len(mean_flux_PX500_39167))), edgecolor='black', color='white')
ax_1_1.axvline(avg_PX500_39167, color='r', linewidth=1, label = '%s Jy'%np.round(avg_PX500_39167, 3))
ax_1_1.margins(x=0)
ax_1_1.set_xlabel(r'$\langle{{E}}\rangle$ (Jy)', fontsize=fontsize)
ax_1_1.legend(loc='center right', fontsize=fontsize, handlelength=1)
ax_1_1.text(0.95, 0.95, 'MJD 58543', transform=ax_1_1.transAxes, fontsize=fontsize, verticalalignment='top', horizontalalignment='right')



plt.savefig("%s_PALL.pdf"%sys.argv[0].split(os.extsep, 1)[0], bbox_inches='tight')
print("%s_PALL.pdf"%sys.argv[0].split(os.extsep, 1)[0])
