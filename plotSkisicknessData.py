# Script to plot skisickness.com-data from
# "Two-Body Kinematics Calculator and Plotter" - 
# http://skisickness.com/2010/04/relativistic-kinematics-calculator/
# For UU course "Accelerators and Detectors, 5c"
#
# PROGRAM DESCRIPTION:
# The program automatically searches for .txt-files in the ./data-folder for data.
# Given it follows generated syntax, it automatically adds it to plottable data.
# Flags can be added to sort the data and keep track of what is what. Labels are
# automatically generated according to filename encoding and flags.
#
# Benjamin Verbeek
# 2022-03-29, Uppsala

# saknar D + 13C och D + 12C och D + ... (alla som har V_c under 2 MeV), D + N? (fler bakgrunds-
# reaktioner). Ligg rent till. Relativ sannolikhet?

# 5 ggr bättre än separationen (5 sigma) (dessutom: VÄLDIGT mkt bakgrund? då kanske 5 sigma för lite... tvärsnitt
# andrav reaktioner?)
# sök specifikation på instrumentet (energiupplösning). Laddad & oladdad i olika!
# neutron mäts på annat sätt...

# skilja på laddningen? Har vi råd? Hur isf, tekniker? 
# hela laddningen ska stoppas upp!
# shift nedåt av stopping i target. 
# instrument ska tåla vakuum.
# energetiskt ger de andra inget.

# Plocka specifikationer vi önskar, sedan: läs

# läs radiation detection and measurement.
# spektroskopisk mätning. Läs de delarna! ej geiger. [Kap 5: jonasitationskammare.]
# testa 3 huvudtekniker. 
# [Fotomultiplikatior / scintillator (kap 9)] funkar med laddade partiklar också! -- BENJAMIN
# [Halvledardetektorer (kap 11) + particle identification (PID)]

# PLOT ONLY CHARGED!!
# Nästa vecka möte: diskutera detta. Bestäm teknik. Sedan: få intro till den tekniken.
# sedan jobb med SP.
# sedan: intro, och vi testar utrustningen. Det jobbet presenteras i rapporten som innehåller:
# - varför tekniken
# - varför specen
# - utvärdering av utrustningen 
# (för både ladd och neutron)

# relativt tvärsnitt ska mätas

# differential cross section: Ungefär samma för de två. 
# Neutronmätning eller inte? p + n eller He3 + T?

# IMPORTS:
import numpy as np
import matplotlib.pyplot as plt
import glob

# filenames whose data you wish to plot:
# find '-E3' for E3, '-E4' for E4, '-th' for theta-corr., '-nc' for no-charge particles
# 'combo' for bothe E3 and E4
# or a specific reaction-code, e.g.: 'dddd' or 'd12cd12c'
display_data = '-nc'
excluded = 'xxx'

MARK = 10   # ticks every n:th data point
# Reactants to look for & how to display them
REACTANTS = {'d':'D', 'n':'n', 'p':'p', 't':'T', '12c':r'^{12}C', '13c':r'^{13}C',\
             '27al':r'^{27}Al', 'h':r'^{3}He', '13n':r'^{13}N', '14n':r'^{14}N', '14c':r'^{14}C'}

if display_data == 'both':
    sort_by = '-E3'
    sort_by2 = '-E4'
else: 
    sort_by = sort_by2 = display_data


# findReactants(encoded string, reactant dictionary):
#   Given the skisickness-encoded filename, extracts reactants for label
#   and returns them in order in display-format (dict value), given a 
#   dictionary of reactants to look for. Roughly a genexp?
#   PRECONDITION: 2 --> 2 particle reactions only.'
#   RETURNS: List of reactants in display-format, in order.
# EXAMPLE:
#   findReactants('d12ct13c') --> ["D", "${^12}C$", "T", "${^13}C$"]
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
files = glob.glob('./data/*.txt')   # all .txt files in ./data
filenames = [f for f in files if ((f.find(sort_by) != -1 or f.find(sort_by2) != -1) and f.find(excluded) == -1)]
print("Found:", filenames)
# what are you plotting?

xaxes = {'-E3':r"$\theta_3$ [deg]", '-E4':r"$\theta_4$ [deg]", '-th':r"$\theta_3$ [deg]", '-nc':r"$\theta$ [deg]", 'both':r"$\theta$ [deg]"}
yaxes = {'-E3':r"$E_3$ [MeV]", '-E4':r"$E_4$ [MeV]", '-th':r"$\theta_4$ [deg]", '-nc':r"$E$ [MeV]", 'both':r"$E$ [MeV]"}

try: xaxis = xaxes[display_data]
except: xaxis = "undefined" # if not covered in dictionary
try: yaxis = yaxes[display_data]
except: yaxis = "undefined" # if not covered in dictionary

i = 0
for filename in filenames:
    # CONSTRUCT LABEL
    start = filename.find('\\p')   # neglect path and leading p
    ending = filename.find("MeV")   # neglect later part
    flags = filename[ending+3:-4]    # fetches flag(s)
    label = filename[start+2:ending].lower()  # all lowercase
    print(filename, label, flags)
    energy = label[-1]  # last number indicates energy in MeV
    r1, r2, r3, r4 = findReactants(label[:-1])
    label = f"${r1} + {r2} \\rightarrow {r3} + {r4}$" + f" [{flags}]"   # with flag

    # READ DATA
    with open(filename) as f:
        lines = f.readlines()

    data = [l.split() for l in lines]
    data = data[2:-2]   # two rows of text, two blank lines.
    data = [[float(v) for v in point] for point in data]    # str to floats
    data = zip(*data)   # converts [(x,y)_i] to [(x_i), (y_i)]

    # PLOT
    linestyles = ['-', '-.', '--', ':']
    plt.plot(*data, marker='|', markevery=MARK, label=label, linestyle=linestyles[i//10])
    plt.legend()
    
    i+=1 

plt.title(f'At $E_k = 2$ MeV. Every {MARK}th point marked.')
plt.xlabel(xaxis)
plt.ylabel(yaxis)
plt.show()