
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
from natsort import index_natsorted 
 ## to be able to sort naturally, so 1,2,14,21 instead of 1,14,2,21
from scipy import stats
from matplotlib import pyplot as plt
# =============================================================================


# ==================Import file================================================
# =============================================================================
project = "210528-515F_806R_qPCR_qqqp_BlackSea_Laura"
csv = 'qPCR_CFX/200715-qPCR_gggp_BS_cDNA_for_Laura.csv'
data = pd.read_csv(csv, delimiter=';')
# =============================================================================

# !!! ==================Variables to set=======================================
# =============================================================================
std_copies = 1.04 # copies/µL of the standard *10^x

# ==================Making a standard curve====================================
# =============================================================================
# Extract standard curve data
stdcurve = data[(data["Content"].str.startswith("Std"))]
# Make empty new column
stdcurve['copies/µL'] = ''
for standard in stdcurve.index:
    power = (stdcurve["Sample"].str.split('^', expand=True)[1].astype(float))
    copies = std_copies * 10^ 
    
    
