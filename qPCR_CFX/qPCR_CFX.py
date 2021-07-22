# =============================================================================
# Title:    qPCR_on_CFX analysis tool
# Version:  1.0
# Author:   Maartje Brouwer
# Goal:     Automating qPCR on CFX analysis
# Date:     210624
# =============================================================================

# Include standard dilution series in PCR1 of a batch
# Include a sample_mix (random mix of some samples from the batch) in all PCRs
# of the batch. Preferable at least in fivefold. used for normalization

# !!! Make sure in the CFX manager you have selected content:
# Std for the standard curve samples and Sample: 10^0, 10^1, 10^2... etc.
# Pos Ctrl for the sample_mix 
# Unkn for samples
# NTC for Non Template Controls

# !!! in CFX manager, click [Export], [Custom Export]
# Uncheck the [Include Run Information Header] checkbox
# Check only the boxes for Content, Sample Name and Cq
# Change Export Format to .xlsx
# Click [Export]
# Open in excel, save as .csv otherwise, Cq values are not correct

# Import needed packages=======================================================
import pandas as pd, numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import math
# =============================================================================


# Import file==================================================================
# =============================================================================
project = "Dina"
PCR1 = '//zeus/mmb/molecular_ecology/mollab_team/Projects/2021/2021_Dina/210721-qpcrDina1.csv'
PCR2 = '//zeus/mmb/molecular_ecology/mollab_team/Projects/2021/2021_Dina/210721-qpcrDina2_trial.csv'
PCR1 = pd.read_csv(PCR1, delimiter=';')
PCR2 = pd.read_csv(PCR2, delimiter=';')

PCRs = [PCR1, PCR2]
# =============================================================================


# !!!Variables to set==========================================================
# =============================================================================
std_copies = 4.06 # copies/µL of the standard *10^x
diluted_PCR = 100

# Making a standard curve======================================================
# =============================================================================
# Extract standard curve data
stdcurve = PCR1[(PCR1["Content"].str.startswith("Std"))]
# From the Sample column, extract the power (copies) of the standards
stdcurve['power'] = (
    stdcurve["Sample"].str.split('^', expand=True)[1].astype(float))

# For each standard, calculate + add copies/µL and log_copies to the dataframe
for standard in stdcurve.index:
    power = stdcurve['power'][standard]
    copies = std_copies * 10 ** power
    log_copies = math.log10(copies)
    stdcurve.loc[standard, 'copies/µL'] = copies
    stdcurve.loc[standard, 'log_copies'] = log_copies

# Linear regression + interpolation for standard curve
slope, yintercept, rv, pv, se = stats.linregress(stdcurve["log_copies"], 
                                                 stdcurve["Cq"])
interp = np.linspace(np.min(stdcurve['log_copies']), 
                            np.max(stdcurve['log_copies']), 
                            num=500)
# calculate efficiency
efficiency = (-1+10**(-1/slope))*100

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
# add efficiency to plot
efficiency = "efficiency = " + str(round(efficiency)) + "%"
ax.text(.7, 0.85, efficiency, 
        size=8, color='purple', 
        transform=ax.transAxes
        )

# make layout fit better
plt.tight_layout()



# =============================================================================
# =============================================================================


# Sample calculations==========================================================
# =============================================================================
for PCR in PCRs:
    data = PCR
    
    
    
    # # Extract sample data
    # samples_raw = data[(data["Content"].str.startswith("Unkn"))]
    # # make str from float Cq values (for next step)
    # samples_raw["Cq"] = samples_raw["Cq"].astype(str)
    # # combine duplicate measurements
    # sample_calculations = samples_raw.groupby("Sample")["Cq"].apply('-'.join)
    # # Convert series to dataframe
    # sample_calculations = sample_calculations.to_frame()

# # Check how many reps each sample has (count the join marks '-' +1)
# for sample in sample_calculations.index:
#     reps = sample_calculations["Cq"][sample].count('-') + 1
#     sample_calculations.loc[sample, "reps"] = reps
# # What is the highest amount of reps
# rep_max = int(sample_calculations["reps"].max())
# # a list for possible seperate Cq values, max 5 repetitions
# Cqs = ["Cq1","Cq2","Cq3","Cq4","Cq5"]
# # Get seperate CQ values, amount depending on number of replicates
# sample_calculations[Cqs[:rep_max]]  = (
#         sample_calculations["Cq"].str.split('-', expand=True).astype(float))

# # remove joined Cq and reps column
# sample_calculations = sample_calculations.drop(['Cq', 'reps'], axis=1)

# # Calculate mean Cq values and stdev
# sample_calculations['mean']= sample_calculations.mean(axis=1)
# sample_calculations['stdev']= sample_calculations.iloc[:, sample_calculations.columns!="mean"].std(axis=1)

# # calculate copies/µL in the DNA extract
# for sample in sample_calculations.index:
#     # calculate from std curve formula (10** because usinjg log-copies)
#     copies = 10**((sample_calculations["mean"][sample] - yintercept) / slope)
#     # add to dataframe, use scientific format, 2 decimal points
#     sample_calculations.loc[sample, "extract copies/µL"] = (
#         "{:.2e}".format(copies))

# # save results
# # save standard curve as .png
# plt.savefig("qPCR_CFX/output/standardcurve_"+ project + ".png")
# # save results in excel
# sample_calculations.to_excel("qPCR_CFX/output/results" + project + ".xlsx")

