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
data = 'molecular_tools/190905-Qubit_LAZ19-1_to_LAZ19-40-End_Point_Results_SYBR.csv'
data = pd.read_csv(data, delimiter=';')
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
fig, ax = plt.subplots(dpi=300)
ax.grid(alpha=0.3)
ax.set_title('standardcurve')
#plot HS curve
HS = ax.scatter(stdcurveHS["ng/ul"], stdcurveHS["End RFU"])
ax.set_xlabel('DNA (ng/µl)')
ax.set_ylabel('end RFU')
ax.plot(
        HS_interp, 
        HSintercept + HSslope * HS_interp,
        label='Linear regression',
        linestyle="--"
        )
# plot BR curve
BR = ax.scatter(stdcurveBR["ng/ul"], stdcurveBR["End RFU"])
ax.set_xlabel('DNA (ng/µl)')
ax.set_ylabel('end RFU')
ax.plot(
        BR_interp, 
        BRintercept + BRslope * BR_interp,
        label='Linear regression',
        linestyle="--"
        )
# =============================================================================


# ==================Calculating concentrations=================================
# Make a new dataframe with the HS Samples + RFU only
DNA_concentrations_HS_raw = data.loc[(data["Sample"].str.startswith("HS_"))]
DNA_concentrations_HS = pd.DataFrame()
DNA_concentrations_HS["Sample"] = DNA_concentrations_HS_raw["Sample"]
DNA_concentrations_HS["RFU"] = DNA_concentrations_HS_raw["End RFU"]
# calculate concentration ((RFU - y-intercept) / slope)
DNA_concentrations_HS["[DNA] ng/µL"] = ''
for sample in DNA_concentrations_HS.index:
    concentration = (
        ((DNA_concentrations_HS['RFU'][sample]) - HSintercept) 
        / HSslope
        ) 
    # add concentration to the dataframe
    DNA_concentrations_HS.at[sample,'[DNA] ng/µL'] = concentration
# Sort samples
DNA_concentrations_HS.sort_values(
    by=['Sample'], 
    inplace=True,
    key=lambda x: np.argsort(index_natsorted(DNA_concentrations_HS["Sample"]))
    )     


# Make a new dataframe with the BR Samples + RFU only
DNA_concentrations_BR_raw = data.loc[(data["Sample"].str.startswith("BR_"))]
DNA_concentrations_BR = pd.DataFrame()
DNA_concentrations_BR["Sample"] = DNA_concentrations_BR_raw["Sample"]
DNA_concentrations_BR["RFU"] = DNA_concentrations_BR_raw["End RFU"]
# calculate concentration ((RFU - y-intercept) / slope)
DNA_concentrations_BR["[DNA] ng/µL"] = ''
for sample in DNA_concentrations_BR.index:
    concentration = (
        ((DNA_concentrations_BR['RFU'][sample]) - BRintercept) 
        / BRslope
        )
    # add concentration to the dataframe
    DNA_concentrations_BR.at[sample,'[DNA] ng/µL'] = concentration
# Sort samples
DNA_concentrations_BR.sort_values(
    by=['Sample'], 
    inplace=True,
    key=lambda x: np.argsort(index_natsorted(DNA_concentrations_BR["Sample"]))
    )    
 
# =============================================================================



                                  

