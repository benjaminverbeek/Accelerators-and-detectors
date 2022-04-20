# Detector calibration script
# Imports ASCII file containing the detector calibration data
# and plots the data. Then fits a Gaussian to the right side of 
# the Amiricium peak of 5.486 MeV (85.2%).
#
# Author: Benjamin Verbeek, Uppsala University, 2022-04-20

from xml.etree.ElementInclude import include
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import glob

# Interesting energy interval
LOW_ENERGY = 4850 # keV
HIGH_ENERGY = 6000 # keV

LOW_COUNT = 0
HIGH_COUNT = 310


files = glob.glob('./calibration-data/*.Spe')   # all .txt files in ./data
print(files)

#filename = "./calibration-data/cali-05-ipe450-100-23-4118.Spe"
filename = files[0]

for filename in files:
    print('='*15,f' Plotting for {filename}', '='*15)
    # Open ASCII file containing the detector calibration data
    with open(filename, 'r') as f:
        # Read the file line by line
        lines = f.readlines()
        # Remove \n characters from the end of each line
        lines = [line.rstrip('\n') for line in lines]

    # Find calibration data:
    cal_index = lines.index("$MCA_CAL:")
    # get calibration coefficients (2nd order polynomial)
    cal_data = lines[cal_index+2].split()[:3]
    cal_data = [float(x) for x in cal_data] # convert to float

    # function converting channel number to energy (in keV)
    def channel2energy(channel):
        return cal_data[0] + cal_data[1]*channel + cal_data[2]*channel**2

    # get channel data:
    ch_index = lines.index("$DATA:")
    n_channels = int(lines[ch_index+1].split()[1])
    ch_data = lines[ch_index+2:ch_index+2+n_channels+1]
    # convert data to integers
    ch_data = [int(x) for x in ch_data]
    ch_data = np.array(ch_data)

    # find energy per channel:
    ch_energy = [channel2energy(x+1) for x in range(len(ch_data))]
    # to np array
    ch_energy = np.array(ch_energy)

    # plot data vs energy:
    plt.figure(figsize=(10,6))
    # set axis width
    plt.xlim(LOW_ENERGY, HIGH_ENERGY)
    plt.ylim(LOW_COUNT, HIGH_COUNT)
    plt.plot(ch_energy, ch_data, '.', label=("Raw counts"))
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.title(f"Detector calibration for {filename}")


    # show plot
    #plt.show()

    # find peak locations: (min distance between peaks is 30 bins)
    peaks, _ = find_peaks(ch_data, height=0.05*max(ch_data), distance=30)
    print(peaks) # indices of peaks
    # skip data 25 bins away or more from the peak. Only go to left side of peak due to contamination.
    peak_ends = [p+25 for p in peaks]

    # mark peaks:
    plt.plot(ch_energy[peaks], ch_data[peaks], 'x', label="Detected peaks")
    # mark end of peaks:
    plt.plot(ch_energy[peak_ends], ch_data[peak_ends], 'x', label="Peak cutoff (arbitrary)")

    # center peak is Amiricium
    # FIT GAUSSIAN TO AMIRICIUM PEAK:
    include_left = 3    # number of bins to include on the left side of the peak
    am_peak = peaks[1] - include_left    # index of the peak minus 3 to get some fall-off too
    am_peak_end = peak_ends[1]

    # Gaussian function:
    def gaus(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    # fit gaussian to data:
    p0 = [100, ch_energy[am_peak], 100]   # initial guess
    popt, pcov = curve_fit(gaus, ch_energy[am_peak:am_peak_end], ch_data[am_peak:am_peak_end], p0)

    print(popt)

    # plot fit:
    gauss_range = ch_energy[am_peak-50:am_peak_end+25]
    plt.plot(gauss_range, [gaus(x,*popt) for x in gauss_range], '-', label=f"Gaussian fit of peak-{include_left}:end_peak for Amiricium \n Found: $\mu$ = {popt[1]:.2f} keV, $\sigma$ = {popt[2]:.2f} keV")
    print(f"Found mean: {popt[1]} +/- {np.sqrt(pcov[1,1])}")
    print(f"Found sigma: {popt[2]} +/- {np.sqrt(pcov[2,2])}")
    plt.grid(True)
    plt.legend()
    # save figure as png:
    plt.savefig(f"{filename[:-4]}.png", dpi=400) # can change dpi later.
    plt.show()

