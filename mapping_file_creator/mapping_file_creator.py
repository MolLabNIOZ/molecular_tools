#!/usr/bin/python3
"""
Automated mapping_file creator
VERSION: Feb2022

If you have different primer sets within 1 sequencing lane
make seperate mapping_files for the different primer sets

Fill in a template. this can be .xlsx or .csv:
Template can be found here: 
//zeus.nioz.nl/mmb/molecular_ecology/mollab_team/Sequencing/ngs_sequencing
If .xlsx, sheet with mapping_file data must be called FILL_IN.
$1 and $2 must be Forward_primer and Reverse_primer

You can add as many columns of metadata as you like
In the .xlsx template, add meta data that is equal for all samples to the ReadMe.
In the FILL_IN only enter metadata that differs between samples
This metadata will be added to the end of the mappingfile
The description column will be the last column of the mapping_file.
Anything after that in the template will be disregarded.
Make the description as descriptive as possible (controls etc).

edit:
    220929 Changed barcodes to 3 or 4 digits, so that we can combine
    Linda's primer set with our primer set. Linda's primers will be numbered 
    9001 and on.
    221115 automated NIOZnumber extraction form excel
    240404 Description column doesnt need to be the last column in the excel
    file.

"""

#### Import needed packages
import pandas as pd       # to be able to work with dataframes
from Bio.Seq import Seq   # to be able to do compl_rev
# !!! Set variables for your mappingfile
file_name = 'NIOZ409_Marie_template.xlsx'
# !!! file_path to folder of mapping_file template (.xlsx or .csv)
folder_path = "//zeus/mmb/molecular_ecology/mollab_team/Sequencing/ngs_sequencing/Mapping_files/"
# Change from windows path to unix path
file_path = folder_path + file_name
 
#### Import needed files
if file_path.endswith('.xlsx') or file_path.endswith('.xlsm'):
    sample_file = pd.ExcelFile(file_path)
if file_path.endswith('.csv'):
    sample_file = pd.read_csv(file_path, delimiter=';')

# NIOZ number to name the mapping_file
if file_path.endswith('.xlsx') or file_path.endswith('.xlsm'):
    try:
        ReadMe = sample_file.parse('ProjectInfo')
    except:
        ReadMe = sample_file.parse('ReadMe')
    NIOZnumber = (
        ReadMe.loc[ReadMe['Project_info'] == 'NIOZ_Number', 'Fill_In'].iloc[0])
if file_path.endswith('.csv'):
    NIOZnumber = 'NIOZ???' #!!! fill in yourself

# Only keep FILL_IN sheet
sample_file = sample_file.parse('FILL_IN')

# Creates a list of all the column headers, removes the header named "Description"
# And adds it to the end of the DF.
col_names = sample_file.columns.tolist()
col_names.remove("Description")
col_names.append("Description")
sample_file = sample_file[col_names]
        
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
                'LinkerPrimerSequence', 
                'Barcode_Forward_Primer']], on='Forward_primer', how='left')
df = pd.merge(
    df, 
    rv_primers[['Reverse_primer', 
                'ReversePrimer', 
                'Barcode_Reverse_Primer']], on='Reverse_primer', how='left')

#### Get complement reverse of reverse primer barcode
# Make empty new column
df ['RevComplReverseBarcodesequence'] = ''
# For every sample get the rev_compl of the reverse primer barcode
for sample in df.index:
    Barcode_Reverse_Primer = Seq(df['Barcode_Reverse_Primer'][sample])
    compl_rev_barcode = Barcode_Reverse_Primer.reverse_complement()
    # add reverse complement to the dataframe
    df.at[sample,'RevComplReverseBarcodesequence'] = compl_rev_barcode

#### Construct BarcodeSequence (barcode_fwd + revcompl_barcode_rev)
# Make empty new column
df ['BarcodeSequence'] = ''
# For every sample make the BarcodeSequence
for sample in df.index:
    BarcodeSequence = ((df['Barcode_Forward_Primer'][sample]) + 
                        (df['RevComplReverseBarcodesequence'][sample]))
    # add BarcodeSequence to the dataframe
    df.at[sample,'BarcodeSequence'] = BarcodeSequence

#### Assemble final mappingfile
# Make an empty dataframe
mf = pd.DataFrame()
# First column: SampleID
mf['#SampleID'] = df['#SampleID']
# Next column: BarcodeSequence
mf['BarcodeSequence'] = df['BarcodeSequence']
# Next column: LinkerPrimerSequence 
mf['LinkerPrimerSequence'] = df['LinkerPrimerSequence']
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
# Remove everything after the description column
mf = mf.loc[:,:'Description']

## Save mapping_file
# Save file as RUNID_mapping_file.txt, tab delimited and without the index
mf.to_csv(folder_path + NIOZnumber + "_mapping_file.txt", sep="\t", index=False)

# Saves the ProjectInfo sheet as a .txt file at the same location as the mapping file
ReadMe.to_csv(folder_path + NIOZnumber + "_README!.txt", sep="\t", index=False)

print(f"Your mappingfile can be found here: \n{folder_path}")