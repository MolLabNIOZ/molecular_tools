"""
Script for the Automation of qubit on CFX analysis

!!! Make sure that your sample names start with HS_ or BR_
!!! Make sure the standards are labelled 'BR 0 ng', 'BR 200 ng', 
    'HS 0 ng' and 'HS 20 ng' etc. according to the standard concentration
    
!!! in the CFX maestro software, go to the [end point] tab. Right mouse click
    on the table where the data can be seen and press "Export to Excel..."

221104 SV:  , delimiter added, changed file source and destination to seperate 
            data folders per project
221107 SV:  outcommented the part where HS_ or BR_ were removed because protocol
            was not doing this correctly
240403 MB:  updated so exported excels can be used

"""

#### Import needed packages
import pandas as pd, numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from natsort import index_natsorted 

#### Import file
file_path = 'C:/Users/mbrouwer/Downloads/53_xc.b.cj_Nora1 -  End Point Results.xlsx'
data = pd.read_excel(file_path, 'FAM')

#### Make 2 empty dataframes
DNA_concentrations_BR = pd.DataFrame(columns = ["Sample"])
DNA_concentrations_HS = pd.DataFrame(columns = ["Sample"])
dataframes = [DNA_concentrations_BR, DNA_concentrations_HS]

#### Making HS and BR standard curves
## Make an empty plot with 2 columns
fig, ax = plt.subplots(ncols=2, dpi=300)

## get curve data for different assays
for i, assay in enumerate(['BR', 'HS']):
    # extract standard curve data
    stdcurve = data[(data["Sample"].str.startswith(assay + " "))]
    # If the specific assay was measured, get curve data
    if not stdcurve.empty:
        stdcurve ["ng/µl"] = ''
        # get the concentrations of the original standards
        for standard in stdcurve.index:
            # extract concentration
            concentration = float(stdcurve["Sample"][standard].split(' ')[1])
            # 2µL was added for the standards, so multiply concentration by 2
            concentration = concentration * 2
            # add concentration to the dataframe
            stdcurve.at[standard, "ng/µl"] = concentration
            
        # Get info about the curve
        slope, intercept, rv, pv, se = stats.linregress(stdcurve["ng/µl"].astype(float),
                                                        stdcurve["End RFU"])
        interp = np.linspace(np.min(stdcurve["ng/µl"]),
                              np.max(stdcurve["ng/µl"]),
                              num=500)        
## Plot both standard curves
        # Put the curve of 1st assay on the left, 2nd on the right
        assay_ax = ax[i]
        # Arange title, axes and grid
        x_axis_max = np.max(stdcurve["ng/µl"]) * 1.25
        assay_ax.set_xticks(np.arange(0,x_axis_max,x_axis_max/5)) # x-grid from 0to25, with intervals=5
        assay_ax.grid(alpha=0.3) # transparancy of grid 
        assay_ax.set_title('standardcurve ' + assay, fontsize=16)
        # Plot the dots
        plot = assay_ax.scatter(stdcurve["ng/µl"],   # x-axis
                          stdcurve["End RFU"], # y-axis
                          c='lightsteelblue')  # color of the dots
        assay_ax.set_xlabel('DNA (ng/µl)')
        assay_ax.set_ylabel('end RFU')
        # Draw standard curve line
        assay_ax.plot(interp,                     # x-axis
                (intercept) + (slope) * interp,   # y-axis
                linestyle='--',                   # style of the line
                c = 'cornflowerblue')             # color of the line
        ## Add equation to plot
        equation = "y = " + str(round(slope)) + "X + " + str(round(intercept))
        assay_ax.text(0.1, 0.9, equation, size=8, color='purple', transform=assay_ax.transAxes)

#### Calculating concentrations
        # Get data per assay
        DNA_concentrations = data.loc[(data["Sample"].str.startswith(assay + "_"))]
        # Add sample names and RFU columns to results dataframe
        dataframes[i]["Sample"] = DNA_concentrations["Sample"]
        dataframes[i][assay + "_RFU"] = DNA_concentrations["End RFU"]
        # Add a new column for the concentration
        dataframes[i][assay + "_[DNA] ng/µL"] = ''
        # for every sample calculate concentration ((RFU - y-intercept) / slope)
        for sample in dataframes[i].index:
            concentration = (
                ((dataframes[i][assay + '_RFU'][sample]) - intercept) 
                / slope) 
            # Add concentration to the dataframe
            dataframes[i].at[sample, assay + '_[DNA] ng/µL'] = concentration
            # Get rid of assay name in sample names, to be able to merge samples
            dataframes[i]['Sample'][sample] = (dataframes[i]['Sample'][sample].split(assay + '_')[1])

## Make plot layout fit better
plt.tight_layout()
## Save standard curve as .png
plt.savefig(file_path.replace(".xlsx","_standardcurve.png"))

# Merge HS and BR measurements
DNA_concentrations = pd.merge(
    DNA_concentrations_HS, 
    DNA_concentrations_BR,
    how = 'outer',      # Use union of keys from both frames
    indicator = 'merge' # new column that tells if both HS and BR are measured
    )
# Sort samples, naturally (so 1,2,14,21 instead of 1,14,2,21)
DNA_concentrations.sort_values(
    by=['Sample'], 
    inplace=True,
    key=lambda x: np.argsort(index_natsorted(DNA_concentrations["Sample"]))
    )

# Final column for concentration, chose if both are measured for BR or HS
DNA_concentrations['[DNA] ng/µL'] = ''
for index in DNA_concentrations.index:
    # If only HS is measured and concentration is below 5
    # concentration is HS measurement
    if DNA_concentrations.loc[index,'merge'] == 'left_only':
        if DNA_concentrations.loc[index,'HS_[DNA] ng/µL'] < 5:
            DNA_concentrations.loc[index,'[DNA] ng/µL'] = (
                DNA_concentrations.loc[index,'HS_[DNA] ng/µL']
                )
    # If concentration is above 5, advise to measure with BR kit
        else:
            DNA_concentrations.loc[index,'[DNA] ng/µL'] = "measure with BR kit"

    # If only BR is measured and concentration is above 5
    # concentration is BR measurement
    elif DNA_concentrations.loc[index,'merge'] == 'right_only':
        if DNA_concentrations.loc[index,'BR_[DNA] ng/µL'] > 5:
            DNA_concentrations.loc[index,'[DNA] ng/µL'] = (
                DNA_concentrations.loc[index,'BR_[DNA] ng/µL']
                )
    # If concentration is below 5, advise to measure with HS kit
        else:
            DNA_concentrations.loc[index,'[DNA] ng/µL'] = "measure with HS kit"
    
    # If both HS and BR is measured and BR[DNA] > 5 and HS[DNA] > 5,
    # concentration is BR measurement
    elif (DNA_concentrations.loc[index,'merge'] == 'both'
          and DNA_concentrations.loc[index,'BR_[DNA] ng/µL'] > 5
          and DNA_concentrations.loc[index,'HS_[DNA] ng/µL'] > 5
          ):
        DNA_concentrations.loc[index,'[DNA] ng/µL'] = (
            DNA_concentrations.loc[index,'BR_[DNA] ng/µL']
            )

    # If both HS and BR is measured and BR[DNA] <= 5,
    # concentration is HS measurement
    elif (DNA_concentrations.loc[index,'merge'] == 'both' 
          and DNA_concentrations.loc[index,'BR_[DNA] ng/µL'] <= 5
          ):
        DNA_concentrations.loc[index,'[DNA] ng/µL'] = (
            DNA_concentrations.loc[index,'HS_[DNA] ng/µL']
            )

    else:
        DNA_concentrations.loc[index,'[DNA] ng/µL'] = "indecisive"

# Remove the merge column
DNA_concentrations.drop(columns=['merge'], inplace=True)

# Save results in excel file
DNA_concentrations.to_excel(file_path.replace(".xlsx","_analyses.xlsx"), index=False)