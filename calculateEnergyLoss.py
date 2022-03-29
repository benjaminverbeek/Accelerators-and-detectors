# Goal:
    # Q4: Determine the energy loss of 2 MeV D ions in a 
    # 5.6 μm thick CD2 target and a 1.9 μm thick natC target. 
    # Tables of stopping powers are available in the folder stopping_powers.

from tracemalloc import start
import numpy as np
import pandas as pd

def load_stopping_power(particle, target):
    # Load file
    filename = f"stopping_powers/{particle}_in_{target}.txt"
    df = pd.read_csv(filename,sep='\t')

    # Extract the first and last row index of the data (the data is surrounded by '-')
    start_idx = 0
    end_idx = 0
    for idx, row in df.iterrows():
        if row.values[0][0] == '-':
            if start_idx == 0:
                start_idx = idx + 1
            else:
                end_idx = idx

    # Extract data from data rows
    rows = np.array([])
    for idx, row in df.iterrows():
        if start_idx <= idx < end_idx:
            rows = np.append(rows, row.values[0])
            
    # Extract energy, dE/dx Elec., dE/dx Nuclear
    vals = np.array([np.array([val.strip() for val in row.split(" ") if val != '']) for row in rows])
    energy_list = np.array([float(row[0]) for row in vals]) # MeV
    dEdx_Elec_list = np.array([float(row[2]) for row in vals]) * 1.0600e-1 # MeV / μm
    dEdx_Nuclear_list = np.array([float(row[3]) for row in vals]) * 1.0600e-1 # MeV / μm

    return energy_list, dEdx_Elec_list, dEdx_Nuclear_list



particles = {
    "D": "2H"
}

targets_and_thickness = {
    "CD2-5.6": ("CD2", 5.6),
    "natC-1.9": ("C", 1.9)
}

initial_energy = 2 # MeV
particle = particles["D"]
target, thickness = targets_and_thickness["CD2-5.6"]
step_size = 0.01 # μm
steps = int(thickness // step_size)

# Load stopping power data
energy_list, dEdx_Elec_list, dEdx_Nuclear_list = load_stopping_power(particle, target)

# Traverse through the medium
energy = initial_energy
traversed_length = 0
for i in range(steps):
    # # Find nearest energy that matches the current energy, and
    # # find the index of minimum element from the differencearray
    # difference_array = np.absolute(energy_list - energy)
    # index = difference_array.argmin()

    # Interpolate the stopping power to the energy
    interp_dEdx_Elec = np.interp(energy, energy_list, dEdx_Elec_list)
    interp_dEdx_Nuclear = np.interp(energy, energy_list, dEdx_Nuclear_list)

    # Sum stopping powers (with a total minus sign)
    total_stopping_power = -(interp_dEdx_Elec + interp_dEdx_Nuclear)

    print(f"Step {i} (traversed {(i*step_size):.2f} μm): energy {energy:.4f} MeV, total stopping power {total_stopping_power:.4f}")
    # Recalculate the energy
    energy = energy + total_stopping_power * step_size

    
print(f"Particle: {particle}, target: {target}, thickness {thickness} μm")
print(f"\tFinal energy is {energy:.4f} MeV. Energy loss is {(initial_energy - energy):.4f} MeV")
