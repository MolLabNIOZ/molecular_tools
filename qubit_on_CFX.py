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
# Make a new dataframe with the HS Samples
DNA_concentrations_HS = data.loc[(data["Sample"].str.startswith("HS_")), "Sample"]
# Add the measured RFU
for sample in DNA_concentrations_HS.Sample:
    RFU = (data['End RFU'][sample])
    # add BarcodeSequence to the dataframe
    DNA_concentrations_HS.at[sample,'End RFU'] = RFU
    

#  DNA_concentrations_BR = data.loc[(data["Sample"].str.startswith("BR_")), "Sample"]

 
# =============================================================================



                                  

