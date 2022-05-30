'''
Program analyzing the data gathered from the tandem D+D fusion
experiment coundcuted for the course Accelerators and Detectors.

The program reads the data semi-hard coded, and then calibrates
the CD2 data to the C background measurement. This produces our
output signal where we identify relevant peaks and extract the
number of counts in each peak (peak bin indices are hard-coded
by trial-and-error.)

After this, plots are generated of the raw data (with relevant 
x- and y-limits) as well as of the analyzed data.

Benjamin Verbeek
Uppsala University
2022-04-26
'''

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import glob
import seaborn as sns
from scipy.optimize import curve_fit



############## LOAD IN DATA ##############
## FETCH CD2 DATA
print('='*15, f' FETCHING DATA ', '='*15)
files = glob.glob('./lab-data/Group B/*.mpa') # all .mpa files in ./lab-data
#print(files)

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
data_A_CD2 = np.array([int(d) for d in data_A_CD2]) # A
data_B_CD2 = data[data_indices[1][0]:data_indices[1][0] + data_indices[1][1]]
data_B_CD2 = np.array([int(d) for d in data_B_CD2]) # B
data_N_CD2 = data[data_indices[2][0]:data_indices[2][0] + data_indices[2][1]]
data_N_CD2 = np.array([int(d) for d in data_N_CD2]) # N
data_BN_coincidence_CD2 = data[data_indices[3][0]:data_indices[3][0] + data_indices[3][1]]
data_BN_coincidence_CD2 = [int(d) for d in data_BN_coincidence_CD2] # BN coincidence
data_BN_coincidence_CD2 = np.array(data_BN_coincidence_CD2).reshape(int(np.sqrt(data_indices[3][1])), int(np.sqrt(data_indices[3][1])))

## FETCH C DATA
C_data = files[2]
print(C_data)
with open(C_data) as f:
    data = f.readlines()
    # Remove \n characters from the end of each line
    data = [line.rstrip('\n') for line in data]

data_indices = [(i+2, int(s[len('[DATA0,'):-2])-1) for i, s in enumerate(data) if '[DATA' in s]
data_indices.append( (data.index('[CDAT0,262144 ]')+1, 262144) )
print(data_indices)

data_A_C = data[data_indices[0][0]:data_indices[0][0] + data_indices[0][1]]
data_A_C = np.array([int(d) for d in data_A_C]) # A

data_B_C = data[data_indices[1][0]:data_indices[1][0] + data_indices[1][1]]
data_B_C = np.array([int(d) for d in data_B_C]) # B

data_N_C = data[data_indices[2][0]:data_indices[2][0] + data_indices[2][1]]
data_N_C = np.array([int(d) for d in data_N_C]) # N

data_BN_coincidence_C = data[data_indices[3][0]:data_indices[3][0] + data_indices[3][1]]
data_BN_coincidence_C = np.array([int(d) for d in data_BN_coincidence_C]) # BN coincidence, C
data_BN_coincidence_C = data_BN_coincidence_C.reshape(int(np.sqrt(data_indices[3][1])), int(np.sqrt(data_indices[3][1])))
print('='*11, f' FINISHED FETCHING DATA ', '='*11, '\n')
############### END FETCH DATA ###############





############### ANALYSIS ###############
print('='*18, f' ANALYSIS ', '='*18)
## COMPENSATIONS
## CALIBRATION FOR DETECTOR B USING CARBON PEAK AT 4500 keV
carbonpeak = (1000, 1090)   # bin numbers for carbon peak in detector B
tritiumpeak = (715, 801)    # bin numbers for tritium peak in signal B
tritium_protonpeak = (503, 635) # bin numbers for proton (T)-peak in signal A
neutronpeak = (574,615)

n_CD2_carbon = sum(data_B_CD2[carbonpeak[0]:carbonpeak[1]])
n_C_carbon = sum(data_B_C[carbonpeak[0]:carbonpeak[1]])

calibration_B = n_CD2_carbon/n_C_carbon
print("Calibration B:", calibration_B)
# Calibrate data from B detectors:
data_B_C = data_B_C*calibration_B
print("Carbon proportion after calibration: ", sum(data_B_C[carbonpeak[0]:carbonpeak[1]])/sum(data_B_CD2[carbonpeak[0]:carbonpeak[1]]))
signal_B = data_B_CD2 - data_B_C

## CALIBRATION FOR DETECTOR A USING TIME AND CURRENT
## CONSTANTS
avrCD2curr = 0.4127 # nA
avrCcurr = 0.716    # nA
liveCD2time = 14.68 # min
liveCtime = 9   # min

# Compensation factors to calibrate C to CD2
currComp = avrCD2curr/avrCcurr
timeComp = liveCD2time/liveCtime
calibration_timeCurrent = currComp*timeComp
print("Calibration time & current:", calibration_timeCurrent)
# Calibrate data from A detectors:
data_A_C = data_A_C*calibration_timeCurrent
signal_A = data_A_CD2 - data_A_C # remove background

# CALIBRATION FOR NEUTRON DETECTOR:
# Calibrate data from N detectors:
data_N_C = data_N_C*calibration_timeCurrent # copy calibration for A
signal_N = data_N_CD2 - data_N_C


## COUNTING
print("\n")
print("COUNTING:")
n_tritium = sum(signal_B[tritiumpeak[0]:tritiumpeak[1]])
print("#Tritium from detector B (T peak):", n_tritium)

nTritiumReact_A = sum(signal_A[503:635])
print("#Tritium from detector A (p peak):", nTritiumReact_A)

nCoincEvents_BN = sum(sum(data_BN_coincidence_CD2[142:155, 120:255]))
print("#signal coincidence events:", nCoincEvents_BN)

nNeutrons_N = sum(signal_N[neutronpeak[0]:neutronpeak[1]])
print("#Neutrons from detector N:", nNeutrons_N)

distB = 19.5 # cm
distA = 10.2 # cm
distN = 55.5 # cm
areaB = 450 # mm^2
areaA = 100 # mm^2
areaN = 65**2*3.1415

#print(areaN/distN**2, areaB/distB**2) # ---> B is bottleneck. N is "catch-all"
## NEUTRON CALIBRATION FACTOR
neutronEfficiencyFactor = 2
neutronDistanceFactor = (distA/distN)**2    # area scales with square of dist
neutronAreaFactor = areaN/areaA
neutronTotalFactor = 1/neutronEfficiencyFactor * neutronDistanceFactor * neutronAreaFactor
nHe3_byCoinc = nCoincEvents_BN*neutronEfficiencyFactor 
nNeutrons_calibrated = nNeutrons_N/neutronTotalFactor
print("#Calibrated coinc. events (only efficiency factor):", nHe3_byCoinc)
print("Neutron factor detector A vs N:", neutronTotalFactor)
print("Neutrons calibrated to detector A:", nNeutrons_calibrated)


# LIST OF CORRECTIONS:
# 1. Strategy 1: Tritium vs He3 (by coincidence)
#    Calibration of B detector: n_C_carbon/n_CD2_carbon to match 4500 keV C12->C13 peak
#    which is used to find counts in T-peak from detector B. 
#    Further, for this strategy, we compare the T-peak to the coincidence measure. Since
#    the He3 from the coincidence is also measured in detector B, and since size-wise it's
#    the bottleneck, we only correct the coincidence by the neutron detector efficiency.
#
# 2. Strategy 2: We compare neutron counts (from He3) vs proton counts (from Tritium).
#    The protons are counted from detector A, and the neutrons from detector N. For A we
#    correct for runtime and current, and then take into account the relative area and 
#    distance of the respective detectors, as well as the efficiency of N.

# Print results: #protons / #neutrons and #Tritium / #He3
print("\n")
print("RESULTS:")
print(f"#protons (T) / #neutrons (He3): {nTritiumReact_A/nNeutrons_calibrated:.3f}")
print(f"#Tritium / #He3: {n_tritium/nHe3_byCoinc:.3f}")


print('='*14, f' FINISHED ANALYSIS ', '='*14, '\n')
############# END ANALYSIS #############





########## PLOTS #############
print('='*15, f' GENERATING PLOTS ', '='*15)
## CD2 TARGET PLOTS
# Plot data_A_CD2 vs index (bin number)
plt.figure(1)
plt.plot([bin2energyA(i) for i in range(data_indices[0][1])], data_A_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector A: CD2 data')
plt.xlim(000,7500)
#plt.ylim(0,20)
plt.savefig('./results/A_CD2.png', format='png')
#plt.yscale('log')
plt.close()

# plot data_B_CD2 vs index (bin number)
plt.figure(2)
plt.plot([bin2energyB(i) for i in range(data_indices[1][1])], data_B_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector B: CD2 data')
plt.xlim(2000,6000)
plt.ylim(0,250)
plt.savefig('./results/B_CD2.png', format='png')
#plt.yscale('log')
plt.close()

# plot data_N_CD2 vs index (bin number)
plt.figure(3)
plt.plot([bin2energyN(i) for i in range(data_indices[2][1])], data_N_CD2)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector N: CD2 data')
plt.xlim(570,630)
plt.ylim(0,3000)
plt.savefig('./results/N_CD2.png', format='png')
#plt.yscale('log')
plt.close()

## COINCIDENCE PLOT
plt.figure(4)
sns.heatmap(data_BN_coincidence_CD2, linewidth=0, norm=LogNorm(), yticklabels=5)
plt.gcf().gca().invert_yaxis()
plt.ylim(80,200)
plt.xlim(0,360)
plt.title('Coincidence between B and N: CD2 data')
plt.savefig('./results/BN_coincidence_CD2.png', format='png')
plt.close()

## CARBON TARGET PLOTS
# Plot data_A_C vs index (bin number)
plt.figure(5)
plt.plot([bin2energyA(i) for i in range(data_indices[0][1])], data_A_C)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector A: C data (compensated for current & runtime)')
plt.xlim(000,7500)
#plt.ylim(0,20)
plt.savefig('./results/A_C.png', format='png')
#plt.yscale('log')
plt.close()

# plot data_B_C vs index (bin number)
plt.figure(6)
plt.plot([bin2energyB(i) for i in range(data_indices[1][1])], data_B_C)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector B: C data (compensated for carbon peak)')
plt.xlim(2000,6000)
plt.ylim(0,250)
plt.savefig('./results/B_C.png', format='png')
#plt.yscale('log')
plt.close()

# plot data_N_C vs index (bin number)
plt.figure(7)
plt.plot([bin2energyN(i) for i in range(data_indices[2][1])], data_N_C)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector N: C data (compensated for current & runtime)')
plt.xlim(570,630)
plt.ylim(0,3000)
plt.savefig('./results/N_C.png', format='png')
#plt.yscale('log')
plt.close()

plt.figure(10)
sns.heatmap(data_BN_coincidence_C, linewidth=0, norm=LogNorm())
plt.gcf().gca().invert_yaxis()
plt.title('Coincidence between B and N: C data')
plt.ylim(80,200)
plt.xlim(0,360)
plt.savefig('./results/BN_coincidence_C.png', format='png')
print("Coinc energies: ",bin2energyB(100*4), bin2energyB(152*4))
plt.close()


## ANALYSIS PLOTS
# normalisera datan med protontoppen från kol
#signal_N = 
# Plot signal_A vs index (bin number)
plt.figure(8)
plt.plot([bin2energyA(i) for i in range(data_indices[0][1])], signal_A)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector A: CD2 minus rescaled C (= signal)')
plt.xlim(000,7500)
#plt.show()
#plt.ylim(0,20)
plt.savefig('./results/signal_A.png', format='png')
#plt.yscale('log')
plt.close()

# plot only data_BN_coincidence_CD2[142:155, 120:255]
plt.figure(11)
sns.heatmap(data_BN_coincidence_CD2[142:155, 120:255], linewidth=0, norm=LogNorm())
plt.gcf().gca().invert_yaxis()
plt.title('Coincidence between B and N: C data, counted')
plt.close()

# plot only signal_A[503:635]
plt.figure(12)
plt.plot(signal_A[503:635])
plt.xlabel('Bin number')
plt.ylabel('Counts')
plt.title('Detector A: Counted signal')
plt.close()
#plt.show()

# Plot signal B
plt.figure(13)
plt.plot([bin2energyB(i) for i in range(data_indices[1][1])], signal_B)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector B: CD2 minus carbon calibrated C (= signal)')
plt.xlim(2100,6000)
plt.ylim(-125,250)
plt.savefig('./results/signal_B.png', format='png')
plt.close()

# plot signal N
plt.figure(14)
plt.plot([bin2energyN(i) for i in range(data_indices[2][1])], signal_N)
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.title('Detector N: CD2 minus carbon calibrated C (= signal)')
plt.xlim(570,630)
plt.ylim(0,3000)
plt.savefig('./results/signal_N.png', format='png')
#plt.show()
plt.close()
########## END PLOTS #############


#################################

# Gaussian function:
def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

lower, upper = 580, 627
fitLow, fitUp = 450, 660

# fit gaussian to signal_A[560:635]
popt, pcov = curve_fit(gauss, np.arange(lower,upper), signal_A[lower:upper], p0=[1, 590, 10])
print("Gaussian fit:", popt)
# plot fit and data
plt.figure(13)
plt.plot(np.arange(fitLow,fitUp), signal_A[fitLow:fitUp], label='data')
plt.plot(np.arange(fitLow,fitUp), gauss(np.arange(fitLow,fitUp), *popt), label='fit')
# plot vertical bars at lower, upper
plt.vlines(lower, 0, gauss(lower,*popt), label='lower', colors='r')
plt.vlines(upper, 0, gauss(upper,*popt), label='upper', colors='r')
plt.legend()
plt.title('Gaussian fit to right hand side of signal_A')
plt.savefig('./results/gauss_fit_sigA.png', format='png')

fittedTcount = sum(gauss(np.arange(fitLow,fitUp), *popt))
print("Fitted T count:", fittedTcount)

print(f"Fitted tritiums = {fittedTcount}, calibrated neutrons = {nCoincEvents_BN*neutronTotalFactor}, ratio = {fittedTcount/(nCoincEvents_BN*neutronTotalFactor)}")

#plt.show()


# TODO: Ta alla neutroner jmf. med protonerna
# TODO: Ta tritiumtoppen och jämför med He3 coincidans.
# TODO: Normalisera med koltoppen från detektor B (ska ha samma antal counts)
# TODO: List all corrections for the two different methods.