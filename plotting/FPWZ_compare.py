import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive
import os
from scipy.signal import find_peaks

# python pol_waterfall.py file Polarisation
# where Polarisation is either I (total intensity), SI (Stokes I), SQ (Stokes Q), SU (Stokes U), L (linear sqrt(SQ^2+SU^2)), SV (Stokes V)

a = psrchive.Archive_load(sys.argv[1])
# polarisation type I,SI,SQ,SU,L,SV
p = sys.argv[2]
# neutron star period (ms)
period = 1/float(sys.argv[3])

### DETERMINE PEAK FLUX AND INDEX FOR PLOT CENTRING
c = a.clone()
c.remove_baseline()
c.tscrunch()
c.fscrunch()
c.pscrunch()
data = c.get_data()
nsub, npol, nchan, nbin = data.shape


# peak and index
flux = data[0,0,0,:]
peaks, _ = find_peaks(flux)
peak_flux = np.sort(flux[peaks])[-1]
peak_idx = np.where(flux==peak_flux)[0][0]

# on-pulse phase start and finish
p1 = np.round(peak_idx/nbin - 0.1, 4)
p2 = np.round(peak_idx/nbin + 0.1, 4)
# off-pulse phase start and finish
p1_off = 0
p2_off = 0.2
# zoom by phase on each side
dp = 0.0985


if sys.argv[2] == "I":
    ### ZOOMED POLARISATION
    c1 = a.clone()
    c1.remove_baseline()
    c1.tscrunch()
    #c1.fscrunch(4)
    c1.pscrunch()
    data1 = c1.get_data()
    data1[0,0,0:43,:] = 0
    nsub, npol, nchan, nbin = data1.shape

    # on-pulse phase bin start and finish
    ps = int(np.round(p1*nbin))
    pf = int(np.round(p2*nbin))
    ps_off = int(np.round(p1_off*nbin))
    pf_off = int(np.round(p2_off*nbin))
    # on-pulse phase bin start and finish
    p3 = p1+dp
    p4 = p2-dp
    p3_off = p1_off+dp
    p4_off = p2_off-dp
    # on-pulse phase start and finish
    psz = int(np.round(p3*nbin))
    pfz = int(np.round(p4*nbin))
    psz_off = int(np.round(p3_off*nbin))
    pfz_off = int(np.round(p4_off*nbin))

    # zoomed intensity
    Iz = data1[0,0,:,psz:pfz]
    Iz_off = data1[0,0,:,psz_off:pfz_off]


    ### POLARISATION
    c1.bscrunch(8)
    data2 = c1.get_data()
    data2[0,0,0:43,:] = 0
    nsub, npol, nchan, nbin = data2.shape

    # on-pulse phase bin start and finish
    ps = int(np.round(p1*nbin))
    pf = int(np.round(p2*nbin))
    ps_off = int(np.round(p1_off*nbin))
    pf_off = int(np.round(p2_off*nbin))

    # intensity
    I = data2[0,0,:,ps:pf]
    I_off = data2[0,0,:,ps_off:pf_off]

else:
    ### ZOOMED POLARISATION
    c1 = a.clone()
    c1.remove_baseline()
    c1.tscrunch()
    #c1.fscrunch(4)
    data1 = c1.get_data()
    data1[0,0,0:43,:] = 0
    nsub, npol, nchan, nbin = data1.shape

    # on-pulse phase bin start and finish
    ps = int(np.round(p1*nbin))
    pf = int(np.round(p2*nbin))
    ps_off = int(np.round(p1_off*nbin))
    pf_off = int(np.round(p2_off*nbin))
    # on-pulse phase bin start and finish
    p3 = p1+dp
    p4 = p2-dp
    p3_off = p1_off+dp
    p4_off = p2_off-dp
    # on-pulse phase start and finish
    psz = int(np.round(p3*nbin))
    pfz = int(np.round(p4*nbin))
    psz_off = int(np.round(p3_off*nbin))
    pfz_off = int(np.round(p4_off*nbin))

    # zoomed polarisations
    SIz = data1[0,0,:,psz:pfz]
    SIz_off = data1[0,0,:,psz_off:pfz_off]
    SQz = data1[0,1,:,psz:pfz]
    SQz_off = data1[0,1,:,psz_off:pfz_off]
    SUz = data1[0,2,:,psz:pfz]
    SUz_off = data1[0,2,:,psz_off:pfz]
    Lz = np.sqrt(data1[0,1,:,psz:pfz]**2+data1[0,2,:,psz:pfz]**2)
    Lz_off = np.sqrt(data1[0,1,:,psz_off:pfz_off]**2+data1[0,2,:,psz_off:pfz_off]**2)
    SVz = data1[0,3,:,psz:pfz]
    SVz_off = data1[0,3,:,psz_off:pfz_off]

    ### POLARISATION
    c1.bscrunch(8)
    data2 = c1.get_data()
    data2[0,0,0:43,:] = 0
    nsub, npol, nchan, nbin = data2.shape

    # on-pulse phase bin start and finish
    ps = int(np.round(p1*nbin))
    pf = int(np.round(p2*nbin))
    ps_off = int(np.round(p1_off*nbin))
    pf_off = int(np.round(p2_off*nbin))

    # polarisations
    SI = data2[0,0,:,ps:pf]
    SI_off = data2[0,0,:,ps_off:pf_off]
    SQ = data1[0,1,:,ps:pf]
    SQ_off = data1[0,1,:,ps_off:pf_off]
    SU = data1[0,2,:,ps:pf]
    SU_off = data1[0,2,:,ps_off:pf_off]
    L = np.sqrt(data2[0,1,:,ps:pf]**2+data2[0,2,:,ps:pf]**2)
    L_off = np.sqrt(data2[0,1,:,ps_off:pf_off]**2+data2[0,2,:,ps_off:pf_off]**2)
    SV = data2[0,3,:,ps:pf]
    SV_off = data2[0,3,:,ps_off:pf_off]


#### PLOTTING ####
plt.figure(figsize=(30,15),dpi=300)

### PLOT ZOOMED POLARISATION
xticks = np.round(np.linspace(p3*period,p4*period,num=11),4)
xticks_x = np.linspace(0,pfz-psz,num=len(xticks))
yticks = np.linspace(704,4032, num=14).astype(int)
yticks_y = np.linspace(0,nchan, len(yticks))

ax3 = plt.subplot(325)
plt.imshow(eval(p+'z'), cmap='Spectral', vmin=np.min(eval(p)), vmax=np.max(eval(p)), aspect='auto', origin='lower')
plt.xticks(xticks_x, xticks)
plt.yticks(yticks_y, yticks)
plt.xlabel('Phase')
plt.ylabel('Frequency (MHz)')

xticks = np.round(np.linspace(p3_off*period,p4_off*period,num=11),4)
xticks_x = np.linspace(0,pfz_off-psz_off,num=len(xticks))

ax3_2 = plt.subplot(326)
plt.imshow(eval(p+'z_off'), cmap='Spectral', vmin=np.min(eval(p)), vmax=np.max(eval(p)), aspect='auto', origin='lower')
plt.xticks(xticks_x, xticks)
plt.yticks(yticks_y, yticks)
plt.xlabel('Phase')
plt.ylabel('Frequency (MHz)')


### PLOT POLARISATION
xticks = np.round(np.linspace(p1*period,p2*period,num=11),4)
xticks_x = np.linspace(0,pf-ps,num=len(xticks))
yticks = np.linspace(704,4032, num=14).astype(int)
yticks_y = np.linspace(0,nchan, len(yticks))

ax2 = plt.subplot(323)
plt.imshow(eval(p), cmap='Spectral', vmin=np.min(eval(p)), vmax=0.3*np.max(eval(p)), aspect='auto', origin='lower')
plt.xticks(xticks_x, xticks)
plt.yticks(yticks_y, yticks)
plt.ylabel('Frequency (MHz)')

xticks = np.round(np.linspace(p1_off*period,p2_off*period,num=11),4)
xticks_x = np.linspace(0,pf_off-ps_off,num=len(xticks))

ax2_2 = plt.subplot(324)
plt.imshow(eval(p+'_off'), cmap='Spectral', vmin=np.min(eval(p)), vmax=0.3*np.max(eval(p)), aspect='auto', origin='lower')
plt.xticks(xticks_x, xticks)
plt.yticks(yticks_y, yticks)
plt.ylabel('Frequency (MHz)')


### PLOT FLUX DENSITY
nsub, npol, nchan, nbin = data.shape

# on-pulse phase start and finish
ps = int(np.round(p1*nbin))
pf = int(np.round(p2*nbin))
ps_off = int(np.round(p1_off*nbin))
pf_off = int(np.round(p2_off*nbin))

xticks = np.round(np.linspace(p1*period,p2*period,num=11),4)
xticks_x = np.linspace(0,pf-ps,num=len(xticks))

ax1 = plt.subplot(321)
plt.plot(data[0,0,0,ps:pf], c='black')
ax1.set_xlim(0,pf-ps)
ax1.set_ylim(-1000,peak_flux+1000)
plt.xticks(xticks_x, xticks)
plt.ylabel('Flux Density (mJy)')
plt.title('On-pulse %s %s'%(p,sys.argv[1]))

xticks = np.round(np.linspace(p1_off*period,p2_off*period,num=11),4)
xticks_x = np.linspace(0,pf_off-ps_off,num=len(xticks))

ax1_2 = plt.subplot(322)
plt.plot(data[0,0,0,ps_off:pf_off], c='black')
ax1_2.set_xlim(0,pf_off-ps_off)
ax1_2.set_ylim(-1000,peak_flux+1000)
plt.xticks(xticks_x, xticks)
plt.ylabel('Flux Density (mJy)')
plt.title('Off-pulse %s %s'%(p,sys.argv[1]))

### SAVE FIGURE
plt.savefig('f_pol_waterfall_comp_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
print('f_pol_waterfall_comp_%s_%s.pdf'%(p,sys.argv[1].split(os.extsep, 1)[0]))
