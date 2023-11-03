"""
This script uses tapestation results for calculating values for equimolar 
pooling.

Input should be a compactRegionTable as produced by the tapestation (.csv) and
the total amount of DNA (ng) you want in your final pool.

Output can either be the recommendation to dilute the samples first, or to
directly pool the undiluted PCR products.
Furthermore, output will be a list with volumes to pool and if necesarry
the volumes of water/DNA needed to dilute the samples.
"""

# Import needed packages=======================================================
import pandas as pd
# =============================================================================

# Variables to set ============================================================
#### Where is the compactRegionTable .csv located?
filepath = "C:/Users/mbrouwer/OneDrive - NIOZ/Documenten/GitHub/molecular_tools/tapestation/test_compactRegionTable_2.csv"

#### How much PCR product is available (µL)
PCR_volume = 35

#### How much DNA do you want to send for sequencing? (ng)
total_ng = 1000 
    # The script multiplies this by 2, to take into account you will loose DNA
    # during clean-up
# =============================================================================

# Data analysis ===============================================================
#### Read the compactRegionTable .csv and put into a dataframe
data = pd.read_csv(filepath, encoding='unicode-escape')

#### Get a list with all concentrations
concentrations = data['Conc. [ng/µl]'].tolist()

#### Check based on the sampleset what samples have sufficient DNA to contribute
def remove_insufficient_samples():
    recursive_loop = False
    for concentration in concentrations:
        number_of_samples = len(concentrations)
        if concentration * PCR_volume < (total_ng * 2) / number_of_samples:
            concentrations.remove(concentration)
            recursive_loop = True

        if recursive_loop:
            remove_insufficient_samples()
        else:
            pass
remove_insufficient_samples()
ng_per_sample = (total_ng * 2) / len(concentrations)

#### 


# =============================================================================
#%%

