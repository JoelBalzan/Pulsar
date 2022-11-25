import numpy as np
import sys
import matplotlib.pyplot as plt
import psrchive


a = psrchive.Archive_load(sys.argv[1])
a.remove_baseline()
a.tscrunch()
a.fscrunch()
data = a.get_data()
nsub, npol, nchan, nbin = data.shape


delta = float(sys.argv[2])
linear = np.sqrt(np.power(data[0,1,:,:], 2.) + np.power(data[0,2,:,:], 2.))
## determine the off pulse phase
I = data[0,0,:,:]
I = np.nanmean(I, axis=0)
num = int(nbin*delta)

std = []
for i in range(nbin):
    if (i+num) < nbin:
        temp = I[i:(i+num)]
        std.append(np.nanstd(temp))
    else:
        temp = np.concatenate((I[i:nbin], I[0:(i+num-nbin)]))
        std.append(np.nanstd(temp))

std = np.array(std)
index = np.argmin(std)
baseline = np.nanmean(data[0,:,:,index:(index+num)], axis=-1)

## calculate Q and U rms
sigma = np.nanstd(data[0,:,:,index:(index+num)], axis=-1)

## correcting for linear pol bias
for i in range(nchan):
    mask = linear[i,:] > sigma[0,i]*1.57
    linear[i,mask] = np.sqrt((linear[i,mask]/sigma[0,i])**2. -1.)*sigma[0,i]

    mask = linear[i,:] <= sigma[0,i]*1.57
    linear[i,mask] = 0.

l_baseline = np.nanmean(linear[:,index:(index+nbin)], axis=-1)
for i in range(nchan):
    linear[i,:] = linear[i,:] - l_baseline[i]

l_sig = np.nanstd(linear[:,index:(index+num)], axis=-1)
print("Linear RMS %f"%np.nanmean(l_sig))

## calculating PA
pa = []
err = []
phase_pa = []
for j in range(nbin):
    if np.fabs(linear[0,j]/l_sig[0]) >=4:
        U = data[0,2,0,j]
        Q = data[0,1,0,j]
        sig_U = sigma[2,0]
        sig_Q = sigma[1,0]
        pa.append((np.arctan2(U,Q)/2.)*180/np.pi)
        phase_pa.append(j)

        part1 = Q/(Q*Q+U*U)
        part2 = -U/(Q*Q+U*U)
        err.append(0.5*(180/np.pi)*np.sqrt(part1*part1*sig_U*sig_U+part2*part2*sig_Q*sig_Q))

pa = np.array(pa)
err = np.array(err)
phase_pa = np.array(phase_pa)

print(pa.shape, err.shape)


x = np.arange(pa.shape[0])
#ticks = np.round(np.arange(p1,p2+0.01, step=0.01), 2)
#ticks_x = np.linspace(0,pf-ps,num=len(ticks))

plt.figure(figsize=(15,10),dpi=300)
plt.plot(x, pa, label='PA')
plt.plot(x, err, label='Error')
#plt.xticks(ticks, ticks)
plt.xlabel('Phase')
plt.ylabel('Polarisation Angle')
plt.legend()

plt.savefig('PA.png')
