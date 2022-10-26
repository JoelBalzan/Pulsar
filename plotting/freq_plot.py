import numpy as np
import subprocess
import glob
import os

# pav -Gdp -k freq1,freq2 --publnc -g out.ps/cps file.ext

### PS FILES
all_cal = glob.glob("*.rescaled")
for filename in all_cal:
    out = filename.split(os.extsep, 1)[0]
    cmd = 'pav -Gdp -g %s_freq.ps/cps %s'%(out,filename)
    print (cmd)
    subprocess.call(cmd, shell=True)

### PS TO PDF
all_ps = glob.glob("*.ps")
print(all_ps)
for filename in all_ps:
    cmd = 'ps2pdf %s'%(filename)
    print(cmd)
    subprocess.call(cmd, shell=True)

### COMBINE PDF FILES AND REMOVE PS FILES
cmd = 'pdftk *.pdf cat output freq.pdf && rm -f *.ps'
print(cmd)
subprocess.call(cmd, shell=True)
