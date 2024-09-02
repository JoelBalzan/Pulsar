import psrchive
import os
import glob
import numpy as np
import sys
from tqdm import tqdm

import matplotlib.pyplot as plt

ext = sys.argv[1]
# Get a list of all the archive files in the directory
files = sorted(glob.glob("*.%s"%(ext)))

# Initialize an empty list to store the profiles
profiles = []

#counter = 0
# Loop through each archive file
for file in tqdm(files):
	# Load the archive file
	a = psrchive.Archive_load(file)
	a.remove_baseline()
	a.fscrunch()
	data = a.get_data()
	data = np.squeeze(data)
	
	
	# Append the profile to the list
	profiles.append(data)

	#counter = counter + 1
	#print("%s/%s"%(counter,len(files))," files completed", end='\r')
	
# Calculate the average profile
average_profile = sum(profiles) / len(files)

# Plot the average profile
### PLOTTING ###
A4x, A4y = 8.27, 11.69
fontsize = 10
fig = plt.figure(figsize=(A4x, A4x), dpi=300)

plt.plot(average_profile, color='black', lw=1)
plt.xlabel('Phase (Bins)')
plt.ylabel('Intensity (mJy)')
plt.title('Average Pulse Profile')
plt.savefig('average_profile.png')