import numpy as np
import glob
import psrchive
import pandas as pd

peak_flux = []

for ar in glob.glob("*.rescaled"):
    a = psrchive.Archive_load(ar)
    a.remove_baseline()
    a.tscrunch()
    a.fscrunch()
    a.pscrunch()
    data = a.get_data()
    flux = np.max(data[0,0,0,:])

    f_row = pd.DataFrame([[a,flux]], columns=['Filename', 'PeakFlux_mJy'])
    peak_flux.append(f_row)

peak_flux = pd.concat(peak_flux, ignore_index=True)
peak_flux.sort_values(by=['PeakFlux_mJy'], ascending=False)

peak_flux.to_csv('peak_flux.csv', index=False)


