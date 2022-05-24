# Author: Benjamin Verbeek, Uppsala University, 2022-04-20

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import glob
import seaborn as sns

avrCD2curr = 0.4127 # nA
avrCcurr = 0.716    # nA

liveCD2time = 14.68 # min
liveCtime = 9   # min

# Compensation factors to calibrate C to CD2
currComp = avrCD2curr/avrCcurr
timeComp = liveCD2time/liveCtime
totalComp = currComp*timeComp
print("Total compensation factor:", totalComp)

#######

files = glob.glob('./lab-data/Group B/*.mpa') # all .mpa files in ./lab-data
print(files)

CD2_data = files[1]
print(CD2_data)
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
bin2energyN = lambda bin: offsets[1] + bin*calfacts[1]  # detector N
bin2energyA = lambda bin: offsets[2] + bin*calfacts[2]  # detector A

data_indices = [(i+2, int(s[len('[DATA0,'):-2])-1) for i, s in enumerate(data) if '[DATA' in s]
data_indices.append( (data.index('[CDAT0,262144 ]')+1, 262144) )
print(data_indices)

# Order found by trial and error
data_A_CD2 = data[data_indices[0][0]:data_indices[0][0] + data_indices[0][1]]
data_A_CD2 = [int(d) for d in data_A_CD2] # A
data_B_CD2 = data[data_indices[1][0]:data_indices[1][0] + data_indices[1][1]]
data_B_CD2 = [int(d) for d in data_B_CD2] # B
data_N_CD2 = data[data_indices[2][0]:data_indices[2][0] + data_indices[2][1]]
data_N_CD2 = [int(d) for d in data_N_CD2] # N
data_BN_coincidence_CD2 = data[data_indices[3][0]:data_indices[3][0] + data_indices[3][1]]
data_BN_coincidence_CD2 = [int(d) for d in data_BN_coincidence_CD2] # BN coincidence
data_BN_coincidence_CD2 = np.array(data_BN_coincidence_CD2).reshape(int(np.sqrt(data_indices[3][1])), int(np.sqrt(data_indices[3][1])))


# Plot data_A_CD2 vs index (bin number)
# NOTE: TODO: INCORRECT CALIBRATION FOR A USED NOW.
plt.figure(1)
plt.plot([bin2energyA(i) for i in range(data_indices[0][1])], data_A_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector A: CD2 data')
plt.xlim(000,7500)
#plt.ylim(0,20)
plt.savefig('./results/A_CD2.pdf', format='pdf')
#plt.yscale('log')

# plot data_B_CD2 vs index (bin number)
plt.figure(2)
plt.plot([bin2energyB(i) for i in range(data_indices[1][1])], data_B_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector B: CD2 data')
plt.xlim(2000,6000)
plt.ylim(0,250)
plt.savefig('./results/B_CD2.pdf', format='pdf')
#plt.yscale('log')

# plot data_N_CD2 vs index (bin number)
plt.figure(3)
plt.plot([bin2energyN(i) for i in range(data_indices[2][1])], data_N_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector N: CD2 data')
#plt.yscale('log')

plt.figure(4)
sns.heatmap(data_BN_coincidence_CD2, linewidth=0, norm=LogNorm(), yticklabels=5)
plt.gcf().gca().invert_yaxis()
plt.ylim(80,200)
plt.xlim(0,360)
plt.title('Coincidence between B and N: CD2 data')
plt.savefig('./results/BN_coincidence_CD2.pdf', format='pdf')



#plt.show()


############ C data #############
C_data = files[2]
print(C_data)
with open(C_data) as f:
    data = f.readlines()
    # Remove \n characters from the end of each line
    data = [line.rstrip('\n') for line in data]

# find calibration parameters
offsets = [float(s[len('caloff='):]) for s in data if 'caloff=' in s]
#print(offsets)
calfacts = [float(s[len('calfact='):]) for s in data if 'calfact=' in s]
#print(calfacts)

#bin2energyB = lambda bin: offsets[0] + bin*calfacts[0]  # detector B
#bin2energyA = lambda bin: offsets[1] + bin*calfacts[1]  # detector A
#bin2energyN = lambda bin: offsets[2] + bin*calfacts[2]  # detector N

data_indices = [(i+2, int(s[len('[DATA0,'):-2])-1) for i, s in enumerate(data) if '[DATA' in s]
data_indices.append( (data.index('[CDAT0,262144 ]')+1, 262144) )
print(data_indices)

data_A_C = data[data_indices[0][0]:data_indices[0][0] + data_indices[0][1]]
data_A_C = [int(d)*totalComp for d in data_A_C] # A
# & Rescale calibration to account for relative beam current

data_B_C = data[data_indices[1][0]:data_indices[1][0] + data_indices[1][1]]
data_B_C = [int(d)*totalComp for d in data_B_C] # B

data_N_C = data[data_indices[2][0]:data_indices[2][0] + data_indices[2][1]]
data_N_C = [int(d)*totalComp for d in data_N_C] # N

data_BN_coincidence_C = data[data_indices[3][0]:data_indices[3][0] + data_indices[3][1]]
data_BN_coincidence_C = [int(d) for d in data_BN_coincidence_C] # BN coincidence, C
data_BN_coincidence_C = np.array(data_BN_coincidence_C).reshape(int(np.sqrt(data_indices[3][1])), int(np.sqrt(data_indices[3][1])))


# Plot data_A_C vs index (bin number)
plt.figure(5)
plt.plot([bin2energyA(i) for i in range(data_indices[0][1])], data_A_C)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector A: C data (compensated for current & runtime)')
plt.xlim(000,7500)
#plt.ylim(0,20)
plt.savefig('./results/A_C.pdf', format='pdf')
#plt.yscale('log')

# plot data_B_C vs index (bin number)
plt.figure(6)
plt.plot([bin2energyB(i) for i in range(data_indices[1][1])], data_B_C)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector B: C data (compensated for current & runtime)')
plt.xlim(2000,6000)
plt.ylim(0,250)
plt.savefig('./results/B_C.pdf', format='pdf')
#plt.yscale('log')

# plot data_N_C vs index (bin number)
plt.figure(7)
plt.plot([bin2energyN(i) for i in range(data_indices[2][1])], data_N_C)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector N: C data (compensated for current & runtime)')
#plt.yscale('log')

plt.figure(10)
sns.heatmap(data_BN_coincidence_C, linewidth=0, norm=LogNorm())
plt.gcf().gca().invert_yaxis()
plt.title('Coincidence between B and N: C data')
plt.ylim(80,200)
plt.xlim(0,360)
plt.savefig('./results/BN_coincidence_C.pdf', format='pdf')
print("Coinc energies: ",bin2energyB(100*4), bin2energyB(152*4))


#plt.show()


# ANALYSIS
print("RUNNING ANALYSIS:")
signal_A = np.array(data_A_CD2) - np.array(data_A_C) # remove background
# Plot signal_A vs index (bin number)
plt.figure(8)
plt.plot([bin2energyA(i) for i in range(data_indices[0][1])], signal_A)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector A: CD2 minus rescaled C (= signal)')
plt.xlim(000,7500)
#plt.show()
#plt.ylim(0,20)
plt.savefig('./results/signal_A.pdf', format='pdf')
#plt.yscale('log')

print("COUNTING:")
nTritiumReact_A = sum(signal_A[503:635])
print("nTritium reacts:", nTritiumReact_A)



nNeutrons_BN = sum(sum(data_BN_coincidence_CD2[120:200, 142:155]))
print("# neutrons:", nNeutrons_BN)

distB = 19.5 # cm
distA = 10.2 # cm
distN = 55.5 # cm
areaB = 450 # mm^2
areaA = 100 # mm^2
areaN = 65**2*3.1415

print(areaN/distN**2, areaB/distB**2) # ---> B is bottleneck. N is "catch-all"


neutronEfficiencyFactor = 2
neutronDistanceFactor = (distB/distA)**2    # area scales with square of dist
neutronAreaFactor = areaA/areaB

neutronTotalFactor = neutronEfficiencyFactor * neutronDistanceFactor * neutronAreaFactor
print("Neutron factor:",neutronTotalFactor)
print("Calibrated neutrons:", nNeutrons_BN*neutronTotalFactor)