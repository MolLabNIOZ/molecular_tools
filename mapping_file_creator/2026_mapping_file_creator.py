#!/usr/bin/python3
"""
Automated mapping_file creator
VERSION: Dec 2025

If you have different primer sets within 1 sequencing lane
make seperate mapping_files for the different primer sets

Fill in the template .xlsx
Template can be found here: 
//zeus.nioz.nl/mmb/molecular_ecology/mollab_team/Sequencing/ngs_sequencing/2026_mappingfile_template.xlsx

"""
#### Change this
file_name = 'NIOZ434_template.xlsx'


#### Import needed packages
import pandas as pd       # to be able to work with dataframes
# from Bio.Seq import Seq   # to be able to do compl_rev

#### Where to find the template
folder_path = "//zeus/mmb/molecular_ecology/mollab_team/Sequencing/ngs_sequencing/Mapping_files/"
file_path = folder_path + file_name

#### Import needed files
sample_file = pd.ExcelFile(file_path)
# Extract README info and NIOZnumber to name the mapping_file
ReadMe = sample_file.parse('ProjectInfo')
NIOZnumber = ( ReadMe.loc[ReadMe['Project_info'] == 'NIOZ_Number', 'Fill_In'].iloc[0])
# Only keep FILL_IN sheet
sample_file = sample_file.parse('FILL_IN')
     
#### Generate sampleIDs
# Create a new empty dataframe
df = pd.DataFrame()
# Insert primers
df['Forward_primer']=(sample_file['Forward_primer'].dropna())
df['Reverse_primer']=(sample_file['Reverse_primer'].dropna())


# Extract primernumbers (last 3-4 characters depending on the primer)
try:
    #checks if last 4 characters are intergers, else the 4th digit is an char
    int(sample_file.iloc[1]['Forward_primer'][-4:])
    df['Forward_primer_number'] = (
            (sample_file['Forward_primer'].str.slice(-4)))
    fw_primer = sample_file.iloc[1]['Forward_primer'][:-4]
except ValueError:
    df['Forward_primer_number'] = (
            (sample_file['Forward_primer'].str.slice(-3)))
    fw_primer = sample_file.iloc[1]['Forward_primer'][:-3]

try:
    int(sample_file.iloc[1]['Reverse_primer'][-4:])
    df['Reverse_primer_number'] = (
            (sample_file['Reverse_primer'].str.slice(-4)))    
    rv_primer = sample_file.iloc[1]['Reverse_primer'][:-4]
except ValueError:
    df['Reverse_primer_number'] = (
            (sample_file['Reverse_primer'].str.slice(-3)))
    rv_primer = sample_file.iloc[1]['Reverse_primer'][:-3]

# Generate SampleIDs
df['#SampleID'] = (
    NIOZnumber + '.' + 
    df ['Forward_primer_number'] + '.' + 
    df ['Reverse_primer_number'])

# Import files with primer
fw_primers = pd.read_excel(
    "molecular_tools/mapping_file_creator/primer_lists.xlsx",
    sheet_name=fw_primer, engine='openpyxl')
rv_primers = pd.read_excel(
    "molecular_tools/mapping_file_creator/primer_lists.xlsx",
    sheet_name=rv_primer, engine='openpyxl')

#### Add primer sequence and barcode sequences from database to samples
df = pd.merge(
    df, 
    fw_primers[['Forward_primer', 
                'ForwardPrimer', 
                'Barcode_Forward_Primer']], on='Forward_primer', how='left')
df = pd.merge(
    df, 
    rv_primers[['Reverse_primer', 
                'ReversePrimer', 
                'Barcode_Reverse_Primer']], on='Reverse_primer', how='left')

# #### Get complement reverse of reverse primer barcode
# # Make empty new column
# df ['RevComplReverseBarcodesequence'] = ''
# # For every sample get the rev_compl of the reverse primer barcode
# for sample in df.index:
#     Barcode_Reverse_Primer = Seq(df['Barcode_Reverse_Primer'][sample])
#     compl_rev_barcode = Barcode_Reverse_Primer.reverse_complement()
#     # add reverse complement to the dataframe
#     df.at[sample,'RevComplReverseBarcodesequence'] = compl_rev_barcode

# #### Construct BarcodeSequence (barcode_fwd + revcompl_barcode_rev)
# # Make empty new column
# df ['BarcodeSequence'] = ''
# # For every sample make the BarcodeSequence
# for sample in df.index:
#     BarcodeSequence = ((df['Barcode_Forward_Primer'][sample]) + 
#                         (df['RevComplReverseBarcodesequence'][sample]))
#     # add BarcodeSequence to the dataframe
#     df.at[sample,'BarcodeSequence'] = BarcodeSequence

#### Assemble final mappingfile
# Make an empty dataframe
mf = pd.DataFrame()
# First column: SampleID
mf['#SampleID'] = df['#SampleID']
# Next column: LinkerPrimerSequence 
mf['ForwardPrimerSequence'] = df['ForwardPrimer']
mf['ReversePrimerSequence'] = df['ReversePrimer']
# Next column: Forward and reverse barcode sequence
mf['Forward_barcode'] = df['Barcode_Forward_Primer']
mf['Reverse_barcode'] = df['Barcode_Reverse_Primer']
# Next columna: primer names
mf['ForwardPrimerName'] = df['Forward_primer']
mf['ReversePrimerName'] = df['Reverse_primer']

#### Add metadata from sample_file
for column in sample_file.columns[2:]:
    mf[column] = sample_file[column]

## Save mapping_file
# Save file as RUNID_mapping_file.txt, tab delimited and without the index
mf.to_csv(folder_path + NIOZnumber + "_mapping_file.txt", sep="\t", index=False)

# Saves the ProjectInfo sheet as a .txt file at the same location as the mapping file
ReadMe.to_csv(folder_path + NIOZnumber + "_README!.txt", sep="\t", index=False)

print(f"Your mappingfile can be found here: \n{folder_path}")