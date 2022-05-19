# Author: Benjamin Verbeek, Uppsala University, 2022-04-20

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import glob

files = glob.glob('./lab-data/Group B/*.mpa') # all .mpa files in ./lab-data
print(files)

CD2_data = files[1]
with open(CD2_data) as f:
    data = f.readlines()
    # Remove \n characters from the end of each line
    data = [line.rstrip('\n') for line in data]


# find calibration parameters
offsets = [float(s[len('caloff='):]) for s in data if 'caloff=' in s]
#print(offsets)
calfacts = [float(s[len('calfact='):]) for s in data if 'calfact=' in s]
#print(calfacts)

bin2energyB = lambda bin: offsets[0] + bin*calfacts[0]  # detector B
bin2energyA = lambda bin: offsets[1] + bin*calfacts[1]  # detector A
bin2energyN = lambda bin: offsets[2] + bin*calfacts[2]  # detector N

data_indices = [(i+2, int(s[len('[DATA0,'):-2])-1) for i, s in enumerate(data) if '[DATA' in s]
print(data_indices)

# Order found by trial and error
data_A_CD2 = data[data_indices[0][0]:data_indices[0][0] + data_indices[0][1]]
data_A_CD2 = [int(d) for d in data_A_CD2] # A
data_B_CD2 = data[data_indices[1][0]:data_indices[1][0] + data_indices[1][1]]
data_B_CD2 = [int(d) for d in data_B_CD2] # B
data_N_CD2 = data[data_indices[2][0]:data_indices[2][0] + data_indices[2][1]]
data_N_CD2 = [int(d) for d in data_N_CD2] # N

# Plot data_A_CD2 vs index (bin number)
plt.figure(1)
plt.plot([bin2energyA(i) for i in range(data_indices[0][1])], data_A_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector A: CD2 data')

# plot data_B_CD2 vs index (bin number)
plt.figure(2)
plt.plot([bin2energyB(i) for i in range(data_indices[1][1])], data_B_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector B: CD2 data')

# plot data_N_CD2 vs index (bin number)
plt.figure(3)
plt.plot([bin2energyN(i) for i in range(data_indices[2][1])], data_N_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector N: CD2 data')

plt.show()


############ C data #############
C_data = files[2]
with open(C_data) as f:
    data = f.readlines()
    # Remove \n characters from the end of each line
    data = [line.rstrip('\n') for line in data]

#print(data[:5])

# find calibration parameters
offsets = [float(s[len('caloff='):]) for s in data if 'caloff=' in s]
#print(offsets)
calfacts = [float(s[len('calfact='):]) for s in data if 'calfact=' in s]
#print(calfacts)

bin2energyB = lambda bin: offsets[0] + bin*calfacts[0]  # detector B
bin2energyA = lambda bin: offsets[1] + bin*calfacts[1]  # detector A
bin2energyN = lambda bin: offsets[2] + bin*calfacts[2]  # detector N

data_indices = [(i+2, int(s[len('[data_A_CD2,'):-2])-1) for i, s in enumerate(data) if '[DATA' in s]
print(data_indices)

data_A_C = data[data_indices[0][0]:data_indices[0][0] + data_indices[0][1]]
data_A_C = [int(d) for d in data_A_C] # A
data_B_C = data[data_indices[1][0]:data_indices[1][0] + data_indices[1][1]]
data_B_C = [int(d) for d in data_B_C] # B
data_N_C = data[data_indices[2][0]:data_indices[2][0] + data_indices[2][1]]
data_N_C = [int(d) for d in data_N_C] # N

# Plot data_A_C vs index (bin number)
plt.figure(1)
plt.plot([bin2energyA(i) for i in range(data_indices[0][1])], data_A_C)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector A: C data')

# plot data_B_C vs index (bin number)
plt.figure(2)
plt.plot([bin2energyB(i) for i in range(data_indices[1][1])], data_B_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector B: C data')

# plot data_N_C vs index (bin number)
plt.figure(3)
plt.plot([bin2energyN(i) for i in range(data_indices[2][1])], data_N_C)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector N: C data')

plt.show()
