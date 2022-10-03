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
    220929. Changed barcodes to 4 digits instead of 3, so that we can combine
    Linda's primer set with our primer set. Linda's primers will be numbered 
    9001 and on.

"""

#### Import needed packages
import pandas as pd       # to be able to work with dataframes
from Bio.Seq import Seq   # to be able to do compl_rev

#### !!! Set variables for your mappingfile
file_name = 'test_user_template_for_mappingfile_creatorpy.xlsx'

# file_path to folder of mapping_file template (.xlsx or .csv)
folder_path = "molecular_tools/mapping_file_creator/"
# Change from windows path to unix path
file_path = folder_path + file_name

# !!! NIOZ number to name the mapping_file
NIOZnumber = 'test_run'	
 
#### Import needed files
if file_path.endswith('.xlsx') or file_path.endswith('.xlsm'):
    sample_file = pd.read_excel(file_path, sheet_name='FILL_IN')
if file_path.endswith('.csv'):
    sample_file = pd.read_csv(file_path, delimiter=';')

#### Generate sampleIDs
# Create a new empty dataframe
df = pd.DataFrame()
# Insert primers
df['Forward_primer']=(sample_file['Forward_primer'])
df['Reverse_primer']=(sample_file['Reverse_primer'])
# Extract primernumbers (last 3 characters)
df['Forward_primer_number'] = (
    (sample_file['Forward_primer'].str.slice(-3)))
df['Reverse_primer_number'] = (
    (sample_file['Reverse_primer'].str.slice(-3)))
# Generate SampleIDs
df['#SampleID'] = (
    NIOZnumber + '.' + 
    df ['Forward_primer_number'] + '.' + 
    df ['Reverse_primer_number'])

#### Get primer and barcode sequence info
# Get primer names without barcode number
fw_primer = sample_file.iloc[1]['Forward_primer'][:-3]
rv_primer = sample_file.iloc[1]['Reverse_primer'][:-3]
# Import files with primer
fw_primers = pd.read_excel(
    "molecular_tools/mapping_file_creator/primer_lists.xlsx",
    sheet_name=fw_primer)
rv_primers = pd.read_excel(
    "molecular_tools/mapping_file_creator/primer_lists.xlsx",
    sheet_name=rv_primer)

#### Add primer sequence and barcode sequences from database to samples
df = pd.merge(df, fw_primers)
df = pd.merge(df, rv_primers)

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
