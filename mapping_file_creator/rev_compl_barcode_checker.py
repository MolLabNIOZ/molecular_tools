#### Import needed packages
import pandas as pd       # to be able to work with dataframes
from Bio.Seq import Seq   # to be able to do compl_rev

#### Import needed files
primers = pd.ExcelFile('primer_lists.xlsx')

fwd_primers = primers.parse('515F_Golay')
rev_primers = primers.parse('951R_Golay')

fwd_primers['RevCompl'] = ''
rev_primers['RevCompl'] = ''

for primer in fwd_primers.index:
    barcode = Seq(fwd_primers['Barcode_Forward_Primer'][primer])
    compl_rev_barcode = barcode.reverse_complement()
    fwd_primers.at[primer, 'RevCompl'] = compl_rev_barcode

for primer in rev_primers.index:
    barcode = Seq(rev_primers['Barcode_Reverse_Primer'][primer])
    compl_rev_barcode = barcode.reverse_complement()
    rev_primers.at[primer, 'RevCompl'] = compl_rev_barcode
    
df = pd.DataFrame()
df['fwd_primer'] = fwd_primers['Forward_primer']
df['fwd_barcode'] = fwd_primers['Barcode_Forward_Primer']
df['revcompl_forward_barcode'] = fwd_primers['RevCompl']
df['rev_primer'] = rev_primers['Reverse_primer']
df['rev_barcode'] = rev_primers['Barcode_Reverse_Primer']
df['revcompl_reverse_barcode'] = rev_primers['RevCompl']

df.to_csv("C:/Users/mbrouwer/OneDrive - NIOZ/Documenten/GitHub/molecular_tools/mapping_file_creator/515F_951R_barcodes.txt", sep='\t', index=False)

