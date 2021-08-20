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
# Make sure for follow-up PCRs the threshold is set to the same value as the 
# PCR including standard dilution series

# !!! in CFX manager, click [Export], [Custom Export]
# Uncheck the [Include Run Information Header] checkbox
# Check only the boxes for Content, Sample Name and Cq
# Change Export Format to .xlsx
# Click [Export]
# Open in excel. Add a column on the right called 'dilution' and add how many 
# times diluted samples are. For undiluted samples enter 1
# Save as .csv otherwise, Cq values are not correct

# Import needed packages=======================================================
import pandas as pd, numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import math
from natsort import index_natsorted 
 ## to be able to sort naturally, so 1,2,14,21 instead of 1,14,2,21
# =============================================================================


# Import file==================================================================
# =============================================================================
project = "DinaJuly21"
PCR1 = '//zeus/mmb/molecular_ecology/mollab_team/Projects/2021/2021_Dina/210721-qpcrDina1.csv'
PCR2 = '//zeus/mmb/molecular_ecology/mollab_team/Projects/2021/2021_Dina/210723-qpcrDina2.csv'
PCR3 = '//zeus/mmb/molecular_ecology/mollab_team/Projects/2021/2021_Dina/210723-qpcrDina3.csv'
PCR4 = '//zeus/mmb/molecular_ecology/mollab_team/Projects/2021/2021_Dina/210723-qpcrDina4.csv'

PCR1 = pd.read_csv(PCR1, delimiter=';')
PCR2 = pd.read_csv(PCR2, delimiter=';')
PCR3 = pd.read_csv(PCR3, delimiter=';')
PCR4 = pd.read_csv(PCR4, delimiter=';')

PCRs = [PCR1, PCR2, PCR3, PCR4]

save_location = '//zeus/mmb/molecular_ecology/mollab_team/Projects/2021/2021_Dina/'
# =============================================================================


# !!!Variables to set==========================================================
# =============================================================================
std_copies = 4.06 # copies/µL of the standard *10^x

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


# Normalize data===============================================================
# =============================================================================

##### Normalize data from subsequent PCRs to PCR with standard dilution series
data = pd.DataFrame()

for PCR in PCRs:
    sample_mix = PCR[(PCR['Content'].str.startswith('Pos Ctrl'))]
      ## Extract sample_mix data
    mean_sample_mix = float(sample_mix['Cq'].mean())
      ## calculate mean Cq of sample_mix 
    if 'Std' in PCR.values:
           normalization = mean_sample_mix
      ## If the specified PCR contains the std dilution series, 
      ## use that sample_mix to normalize
    
    normalization_factor = mean_sample_mix - normalization
      ## Calculate normalization_factor for specified PCR
    
    for sample in PCR.index:
        Cq = PCR['Cq'][sample] + normalization_factor
        PCR.loc[sample, 'corrected Cq'] = Cq
          ## Normalize Cq values based on normalization_factor
    
    data = data.append(PCR)

# Sample calculations==========================================================
# =============================================================================    

##### Extract sample data
samples_raw = data[(data["Content"].str.startswith("Unkn"))]
# make str from float Cq values (for next step)
samples_raw["corrected Cq"] = samples_raw["corrected Cq"].astype(str)
# combine duplicate measurements
sample_calculations = (
    samples_raw.groupby("Sample")["corrected Cq"].apply('-'.join))
# Convert series to dataframe
sample_calculations = sample_calculations.to_frame()

# Check how many reps each sample has (count the join marks '-' +1)
for sample in sample_calculations.index:
    reps = sample_calculations["corrected Cq"][sample].count('-') + 1
    sample_calculations.loc[sample, "reps"] = reps
# What is the highest amount of reps
rep_max = int(sample_calculations["reps"].max())
# a list for possible seperate Cq values, max 5 repetitions
Cqs = (["corrected Cq1",
        "corrected Cq2",
        "corrected Cq3",
        "corrected Cq4",
        "corrected Cq5"])
# Get seperate CQ values, amount depending on number of replicates
sample_calculations[Cqs[:rep_max]]  = (
        sample_calculations["corrected Cq"].
        str.split('-', expand=True).astype(float))

# remove joined Cq and reps column
sample_calculations = sample_calculations.drop(['corrected Cq', 'reps'],axis=1)

# Calculate mean Cq values and stdev
sample_calculations['mean']= sample_calculations.mean(axis=1)
sample_calculations['stdev']= (
    sample_calculations.iloc[:, sample_calculations.columns!="mean"].
    std(axis=1))

# Add dilution factor to the dataframe
sample_calculations = pd.merge(
    sample_calculations, data[['Sample', 'dilution']], how='left', on='Sample')
sample_calculations = sample_calculations.drop_duplicates()


# calculate copies/µL in the DNA extract
for sample in sample_calculations.index:
    # calculate from std curve formula (10** because using log-copies)
    copies = 10**((sample_calculations["mean"][sample] - yintercept) / slope)
    # multiply by dilution factor
    # copies = copies * sample_calculations["dilution"][sample]
    # add to dataframe, use scientific format, 2 decimal points
    sample_calculations.loc[sample, "extract copies/µL"] = (
        "{:.2e}".format(copies))
    
# Sort samples, naturally (so 1,2,14,21 instead of 1,14,2,21)
sample_calculations.sort_values(
    by=['Sample'], 
    inplace=True,
    key=lambda x: np.argsort(index_natsorted(sample_calculations["Sample"]))
    )


# save results
# save standard curve as .png
plt.savefig(save_location + "standardcurve_"+ project + ".png")
# save results in excel
sample_calculations.to_excel(save_location + "results" + project + ".xlsx")
