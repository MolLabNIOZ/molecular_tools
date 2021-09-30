# =============================================================================
# Title:    Reverse Complement
# Version:  1.0
# Author:   Maartje Brouwer
# Goal:     Make reverse complement of a column
# Date:     210928
# =============================================================================

# IMPORT STATEMENTS============================================================
# =============================================================================
import pandas as pd       # to be able to work with tables
import numpy as np
from Bio.Seq import Seq   # to be able to do compl_rev
# =============================================================================

# VARIABLES TO SET#!!!=========================================================
# =============================================================================
file_path = 'molecular_tools/reverse_complement/'
sample_file = file_path + '12S_barcodes_minDistances1_20210916.xlsx'

file = pd.read_excel(sample_file)

### New column with reverse complement barcode
file['reverse_complement']=''
file['match_name']=''
for primer in file.index:
    barcode = Seq(file['barcode'][primer])
    reverse_complement = barcode.reverse_complement()
    file.at[primer,'reverse_complement'] = str(reverse_complement)
    
    ###Check if barcode matches any RC barcode
    match = file.index[file['reverse_complement'] == barcode].tolist()
    if len(match) > 0:
        file.loc[primer, 'match_name'] = file['Name'][match[0]]

  
    ###Check if RC barcode matches any barcode
    match = file.index[file['barcode'] == reverse_complement].tolist()
    if len(match) > 0:
        file.loc[primer, 'match_name'] = file['Name'][match[0]]

    
    

    

file.to_csv(file_path + '12S_reverse_complement_barcodes.csv')
    
    
