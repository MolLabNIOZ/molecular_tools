# =============================================================================
# Title:    Qubit_on_CFX analysis tool
# Version:  1.0
# Author:   Maartje Brouwer
# Goal:     Automating qubit on CFX analysis
# Date:     210423
# =============================================================================

# !!! Make sure that your sample names start with HS_ or BR_
# !!! Make sure the standards are labelled 
# 'BR 0 ng', 'BR 200 ng', 'HS 0 ng' and 'HS 20 ng'

# !!! in CFX manager, save the data as excel
# in excel, save as .csv

# ==================Import needed packages=====================================
import pandas as pd, numpy as np
from natsort import index_natsorted 
 ## to be able to sort naturally, so 1,2,14,21 instead of 1,14,2,21
from scipy import stats
from matplotlib import pyplot as plt
# =============================================================================


# ==================Import file================================================
project = "190905_QubitLAZ19-1_to_LAZ19-40"
csv = '190905-Qubit_LAZ19-1_to_LAZ19-40-End_Point_Results_SYBR.csv'
data = pd.read_csv(csv, delimiter=';')
# =============================================================================


# ==================Making a standard curve====================================
# Extract standard curve data
stdcurveHSraw = data[(data["Sample"].str.startswith("HS "))]
stdcurveBRraw = data[(data["Sample"].str.startswith("BR "))]
# Make empty dataframe
stdcurveHS = pd.DataFrame()
stdcurveBR = pd.DataFrame()
# Cut ng/ul from standard name
stdcurveHS["ng/ul"] = stdcurveHSraw["Sample"].str.split(' ', expand=True)[1].astype(float)
stdcurveBR["ng/ul"] = stdcurveBRraw["Sample"].str.split(' ', expand=True)[1].astype(float)
# Get "End RFU" of the standards
stdcurveHS["End RFU"] = stdcurveHSraw["End RFU"]
stdcurveBR["End RFU"] = stdcurveBRraw["End RFU"]

# Linear regression
HSslope, HSintercept, rv, pv, se = stats.linregress(stdcurveHS["ng/ul"], stdcurveHS["End RFU"])
HS_interp = np.linspace(np.min(stdcurveHS["ng/ul"]), np.max(stdcurveHS["ng/ul"]), num=500)

BRslope, BRintercept, rv, pv, se = stats.linregress(stdcurveBR["ng/ul"], stdcurveBR["End RFU"])
BR_interp = np.linspace(np.min(stdcurveBR["ng/ul"]), np.max(stdcurveBR["ng/ul"]), num=500)

# Plot standard curves
# empty plot
fig, ax = plt.subplots(ncols=2, dpi=300)
HSax = ax[0]
HSax.grid(alpha=0.5)
HSax.set_title('standardcurve HS', fontsize=16)
#plot HS curve
HS = HSax.scatter(stdcurveHS["ng/ul"], stdcurveHS["End RFU"], c='lightsteelblue')
HSax.set_xlabel('DNA (ng/µl)')
HSax.set_ylabel('end RFU')
HSax.plot(
        HS_interp, 
        HSintercept + HSslope * HS_interp,
        label='Linear regression',
        linestyle='--',
        c='cornflowerblue'
        )
# plot BR curve
BRax = ax[1]
BRax.grid(alpha=0.5)
BRax.set_title('standardcurve BR', fontsize=16)
BR = BRax.scatter(stdcurveBR["ng/ul"], stdcurveBR["End RFU"], c='cornflowerblue')
BRax.set_xlabel('DNA (ng/µl)')
BRax.set_ylabel('end RFU')
BRax.plot(
        BR_interp, 
        BRintercept + BRslope * BR_interp,
        label='Linear regression',
        linestyle="--",
        c='blue'
        )
# make layout fit better
plt.tight_layout()
# save standard curve
plt.savefig("output/standardcurve"+ project + ".png")

# =============================================================================


# ==================Calculating concentrations=================================
# Make a new dataframe with the HS Samples + RFU only
DNA_concentrations_HS_raw = data.loc[(data["Sample"].str.startswith("HS_"))]
DNA_concentrations_HS = pd.DataFrame()
DNA_concentrations_HS["Sample"] = DNA_concentrations_HS_raw["Sample"]
DNA_concentrations_HS["HS_RFU"] = DNA_concentrations_HS_raw["End RFU"]
# calculate concentration ((RFU - y-intercept) / slope)
DNA_concentrations_HS["HS_[DNA] ng/µL"] = ''
for sample in DNA_concentrations_HS.index:
    concentration = (
        ((DNA_concentrations_HS['HS_RFU'][sample]) - HSintercept) 
        / HSslope
        ) 
    # add concentration to the dataframe
    DNA_concentrations_HS.at[sample,'HS_[DNA] ng/µL'] = concentration
# Sort samples natrually
DNA_concentrations_HS.sort_values(
    by=['Sample'], 
    inplace=True,
    key=lambda x: np.argsort(index_natsorted(DNA_concentrations_HS["Sample"]))
    )     

# Make a new dataframe with the BR Samples + RFU only
DNA_concentrations_BR_raw = data.loc[(data["Sample"].str.startswith("BR_"))]
DNA_concentrations_BR = pd.DataFrame()
DNA_concentrations_BR["Sample"] = DNA_concentrations_BR_raw["Sample"]
DNA_concentrations_BR["BR_RFU"] = DNA_concentrations_BR_raw["End RFU"]
# calculate concentration ((RFU - y-intercept) / slope)
DNA_concentrations_BR["BR_[DNA] ng/µL"] = ''
for sample in DNA_concentrations_BR.index:
    concentration = (
        ((DNA_concentrations_BR['BR_RFU'][sample]) - BRintercept) 
        / BRslope
        )
    # add concentration to the dataframe
    DNA_concentrations_BR.at[sample,'BR_[DNA] ng/µL'] = concentration
# Sort samples naturally
DNA_concentrations_BR.sort_values(
    by=['Sample'], 
    inplace=True,
    key=lambda x: np.argsort(index_natsorted(DNA_concentrations_BR["Sample"]))
    )    

# Get rid of HS / BR in sample name, to be able to merge
DNA_concentrations_HS['Sample'] = (
    (DNA_concentrations_HS['Sample'].str.split('_', expand=True))[1].astype(str))
DNA_concentrations_BR['Sample'] = (
    (DNA_concentrations_BR['Sample'].str.split('_', expand=True))[1].astype(str))

# Merge HS and BR
DNA_concentrations = pd.merge(
    DNA_concentrations_HS, 
    DNA_concentrations_BR,
    how = 'outer',      # Use union of keys from both frames
    indicator = 'merge' # new column, that tells if both HS and BR are measured
    )
# Sort samples, naturally
DNA_concentrations.sort_values(
    by=['Sample'], 
    inplace=True,
    key=lambda x: np.argsort(index_natsorted(DNA_concentrations["Sample"]))
    )

# Final column for concentration chosen if both are measured for BR or HS
DNA_concentrations['[DNA] ng/µL'] = ''
for index in DNA_concentrations.index:
    if DNA_concentrations.loc[index,'merge'] == 'left_only':
        DNA_concentrations.loc[index,'[DNA] ng/µL'] = DNA_concentrations.loc[index,'HS_[DNA] ng/µL']
    elif DNA_concentrations.loc[index,'merge'] == 'right_only':
        DNA_concentrations.loc[index,'[DNA] ng/µL'] = DNA_concentrations.loc[index,'BR_[DNA] ng/µL']
    elif (DNA_concentrations.loc[index,'merge'] == 'both' 
          and DNA_concentrations.loc[index,'BR_[DNA] ng/µL'] <= 5
          ):
        DNA_concentrations.loc[index,'[DNA] ng/µL'] = DNA_concentrations.loc[index,'HS_[DNA] ng/µL']
    elif (DNA_concentrations.loc[index,'merge'] == 'both'
          and DNA_concentrations.loc[index,'BR_[DNA] ng/µL'] > 5
          or DNA_concentrations.loc[index,'HS_[DNA] ng/µL'] > 10
          ):
        DNA_concentrations.loc[index,'[DNA] ng/µL'] = DNA_concentrations.loc[index,'BR_[DNA] ng/µL']
          
          

# Save results
DNA_concentrations.to_excel('output/results' + project + '.xlsx', index=False)
    
# =============================================================================
                             