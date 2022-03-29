# Script to plot skisickness.com-data from
# "Two-Body Kinematics Calculator and Plotter" - 
# http://skisickness.com/2010/04/relativistic-kinematics-calculator/
# For UU course "Accelerators and Detectors, 5c"

# Benjamin Verbeek 2022-03-29, Uppsala.

# IMPORTS:
import numpy as np
import matplotlib.pyplot as plt
import glob

files = glob.glob('./data/*.txt')   # all .txt files in ./data

# Reactants to look for & how to display them
REACTANTS = {'d':'D', 'n':'n', 'p':'p', '12c':r'^{12}C', '13c':r'^{13}C',\
             '27al':r'^{27}Al', 'h':r'^{3}He', '13n':r'^{13}N'}

# findReactants:
#   Given a string and a list of possible reactants (as dict),
#   identifies and returns a list of reactants in order.
#   PRECONDITION: 2 --> 2 particle reactions only.
def findReactants(s, reac=REACTANTS):
    out = []
    for _ in range(4):
        for key, val in reac.items():
            if s.find(key) == 0:
                out.append(val)
                s = s[len(key):]
                break
    if s != "": print("WARNING: Unknown reactant/bad format. Incorrect output.")
    return out
##############

# filenames whose data you wish to plot:
# find '-E3' for E3, '-E4' for E4, '-th' for theta-corr.
filenames = [f for f in files if f.find('-E3') != -1]
print("Found:", filenames)
# what are you plotting?
xaxis = r"$\theta_3$ [deg]"
yaxis = r"$E_3$ [MeV]"


for filename in filenames:
    #filename = "pdddd2MeV3.txt"
    # construct label:
    start = filename.find('\\p')   # neglect path and leading p
    ending = filename.find("MeV")   # neglect later part
    label = filename[start+2:ending].lower()  # all lowercase
    print(filename, label)
    energy = label[-1]  # last number indicates energy in MeV
    r1, r2, r3, r4 = findReactants(label[:-1])
    label = f"${r1} + {r2} \\rightarrow {r3} + {r4}$"

    with open(filename) as f:
        lines = f.readlines()

    data = [l.split() for l in lines]
    data = data[2:-2]   # two rows of text, two blank lines.
    data = [[float(v) for v in point] for point in data]    # str to floats
    data = zip(*data)   # converts [(x,y)_i] to [(x_i), (y_i)]

    plt.plot(*data, marker='|', markevery=5, label=label)
    plt.title('At $E_k = 2$ MeV. Every 5th point marked. ')
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.legend()


plt.show()