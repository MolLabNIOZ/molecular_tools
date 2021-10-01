# =============================================================================
# Title:    Reverse Complement match check
# Version:  1.0
# Author:   Maartje Brouwer
# Goal:     Check reverse complement matches with barcodes
# Date:     210928
# =============================================================================
# Make sure your file is an .xlsx file
# Make sure the column with barcodes is named 'barcodes'
# =============================================================================

# IMPORT STATEMENTS============================================================
# =============================================================================
import pandas as pd       # to be able to work with tables
from Bio.Seq import Seq   # to be able to do compl_rev
# =============================================================================

# VARIABLES TO SET#!!!=========================================================
# =============================================================================
file_path = 'molecular_tools/reverse_complement/'
sample_file = file_path + '12S_barcodes_minDistances1_20210916.xlsx'
# =============================================================================

# IMPORTING DATA===============================================================
# =============================================================================
file = pd.read_excel(sample_file)
# =============================================================================

# ADD COMPLEMENT REVERSE BARCODES AND CHECK FOR MATCHES========================
# =============================================================================
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
# =============================================================================

# SAVE NEW DATAFRAME===========================================================
# =============================================================================
file.to_excel(file_path + '12S_reverse_complement_barcodes.xlsx')
# =============================================================================
    
    
