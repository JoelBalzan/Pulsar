import plotly.express as px
import sys
import os
import numpy as np
import psrchive

filename, file_extension = os.path.splitext(sys.argv[1])

if file_extension == '.ar' or file_extension == '.calib' or file_extension == '.rescaled':
    a = psrchive.Archive_load(sys.argv[1])
    a.remove_baseline()
    a.tscrunch()
    a.pscrunch()
    a = a.get_data()
    a = np.squeeze(a)
    print(np.shape(a))
    #np.save(sys.argv[1] + '.npy', a)

if file_extension == '.npy':
    a = np.load(sys.argv[1])

fig = px.imshow(a, aspect='auto', origin='lower')

fig.show()