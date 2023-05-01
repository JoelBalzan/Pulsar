import numpy as np
import glob
import matplotlib.pyplot as plt
import psrchive


if glob.glob("peak_flux.npy"):
    peak_flux = np.load("peak_flux.npy")
else:
    peak_flux = []
    for ar in glob.glob("*.rescaled"):
        a = psrchive.Archive_load(ar)
        a.remove_baseline()
        a.tscrunch()
        a.fscrunch()
        a.pscrunch()
        data = a.get_data()
        pk_flux = np.max(data[0,0,0,:])/1000
        peak_flux.append(pk_flux)
    peak_flux = np.array(peak_flux)
    np.save("peak_flux.npy", peak_flux)

### PLOTTING ###
A4x, A4y = 8.27, 11.69
fig = plt.figure(figsize=(A4y, A4x), dpi=300)

plt.plot(peak_flux, c='k', marker='o', lw=1)
plt.xlabel("Pulse Number")
plt.ylabel("Peak Flux Density (Jy)")

plt.savefig("peak_flux.pdf", bbox_inches='tight', dpi=300)