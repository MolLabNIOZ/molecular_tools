
# =============================================================================
# Title:    Qubit_on_CFX analysis tool
# Version:  1.0
# Author:   Maartje Brouwer
# Goal:     Automating qubit on CFX analysis
# Date:     210423
# =============================================================================

# !!! Make sure in the CFX manager you have selected:
# Standard for the standard curve samples
#   And name standards 10^0, 10^1, 10^2... etc.
# unknown for samples
# NTC for Non Template Controls

# !!! in CFX manager, click [Export], [Custom Export]
# Uncheck the [Include Run Information Header] checkbox
# Check only the boxes for Content, Sample Name and Cq
# Change Export Format to .xlsx
# Click [Export]
# in excel, save as .csv otherwise, Cq values are not correct

# ==================Import needed packages=====================================
import pandas as pd, numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import math
# =============================================================================


# ==================Import file================================================
# =============================================================================
project = "Tim_qPCR"
csv = 'qPCR_CFX/compare.csv'
data = pd.read_csv(csv, delimiter=';')
# =============================================================================

# !!! ==================Variables to set=======================================
# =============================================================================
std_copies = 4.06 # copies/µL of the standard *10^x

# ==================Making a standard curve====================================
# =============================================================================
# Extract standard curve data
stdcurve = data[(data["Content"].str.startswith("Std"))]
# From the Sample column, extract the power (copies) of the standards
stdcurve['power'] = (
    stdcurve["Sample"].str.split('^', expand=True)[1].astype(float))
# For each standard, calculate + add copies/µL and log_copies to the dataframe
for standard in stdcurve.index:
    power = stdcurve['power'][standard]
    copies = std_copies * 10 ** power
    log_copies = math.log10(copies)
    stdcurve.at[standard, 'copies/µL'] = copies
    stdcurve.at[standard, 'log_copies'] = log_copies
# Linear regression + interpolation for standard curve
slope, yintercept, rv, pv, se = stats.linregress(stdcurve["log_copies"], 
                                                 stdcurve["Cq"])
interp = np.linspace(np.min(stdcurve['log_copies']), 
                            np.max(stdcurve['log_copies']), 
                            num=500)

# plot standard curve
fig, ax = plt.subplots(dpi=300)     # empty plot
ax.set_xticks(np.arange(0,8,1))     # x-grid from 0to8, with intervals=1
ax.grid(alpha=0.3)                  # transparancy of grid
ax.set_title('standardcurve', fontsize=16)
ax.scatter(stdcurve["log_copies"],  # x-axis
           stdcurve["Cq"],          # y-axis
           c='blue')                # color of the dots
ax.set_xlabel('log10 copies')       # x-axis label
ax.set_ylabel('Cq')                 # y-axis label
# plot linear regression
ax.plot(interp,                     # x-axis
        yintercept + slope * interp,# y-axis
        linestyle='--',             # style of the line
        c='cornflowerblue'          # color of the line
        )
# add equation to plot
equation = "y = " + str(round(slope, 3)) + "X + " + str(round(yintercept,2))
ax.text(.7, 0.9, equation, 
        size=8, color='purple', 
        transform=ax.transAxes
        )
# make layout fit better
plt.tight_layout()
# save standard curve as .png
plt.savefig("qPCR_CFX/output/standardcurve_"+ project + ".png")

# ==================Sample calculations========================================
# =============================================================================
# Extract sample data
samples_raw = data[(data["Content"].str.startswith("Unkn"))]
# make str from float Cq values (for next step)
samples_raw["Cq"] = samples_raw["Cq"].astype(str)
# combine duplicate measurements
sample_calculations = samples_raw.groupby("Sample")["Cq"].apply('-'.join)
# Convert series to dataframe
sample_calculations = sample_calculations.to_frame()

# Get seperate CQ values
sample_calculations = (
    sample_calculations["Cq"].str.split('-', expand=True).astype(float))

# Calculate mean Cq values, stdev and variance
    # variance = how far from the mean the individual observations are
    # stdev = square root of the variance, the amount of variation or dispersion


   




