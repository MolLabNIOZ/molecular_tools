#!/usr/bin/python3
"""
Automated mapping_file creator
VERSION: Dec 2025
Also handles M13 primers

If you have different primer sets within 1 sequencing lane
make seperate mapping_files for the different primer sets

Fill in the template .xlsx
Template can be found here: 
//zeus.nioz.nl/mmb/molecular_ecology/mollab_team/Sequencing/ngs_sequencing/2026_mappingfile_template.xlsx

"""
#### Change this
file_name = 'NIOZ431_M13_2025_mappingfile_template.xlsx'

#### Import needed packages
import pandas as pd       # to be able to work with dataframes
import numpy as np
# from Bio.Seq import Seq   # to be able to do compl_rev

#### Where to find the template
folder_path = "//zeus/mmb/molecular_ecology/mollab_team/Sequencing/ngs_sequencing/Mapping_files/"
file_path = folder_path + file_name

#### Import needed files
sample_file = pd.ExcelFile(file_path)
# Extract README info and NIOZnumber to name the mapping_file
ReadMe = sample_file.parse('ProjectInfo')
NIOZnumber = ReadMe.loc[ReadMe['Project_info'] == 'NIOZ_Number', 'Fill_In'].iloc[0]
# Only keep FILL_IN sheet
sample_file = sample_file.parse('FILL_IN')
     
#### Generate sampleIDs
# Create a new empty dataframe
df = pd.DataFrame()
# Insert primers
df['Forward_primer']=(sample_file['Forward_primer'].dropna())
df['Reverse_primer']=(sample_file['Reverse_primer'].dropna())

primer_sequence_library = pd.read_excel("molecular_tools/mapping_file_creator/primer_lists.xlsx", engine='openpyxl', sheet_name=None)

#### Extract primernumbers (last 3-4 characters depending on the primer) and primer names
M13fwd = False
M13rev = False
for sample in df.index:
    Forward = df['Forward_primer'][sample]
    Reverse = df['Reverse_primer'][sample]
    # Extract barcodenumbers
    df.at[sample,'Forward_primer_number'] = Forward_primer_number = ''.join([char for char in ((df['Forward_primer'][sample])[-4:])[::-1] if char.isdigit()])[::-1]
    df.at[sample,'Reverse_primer_number'] = Reverse_primer_number = ''.join([char for char in ((df['Reverse_primer'][sample])[-4:])[::-1] if char.isdigit()])[::-1]
    # Extract primer names    
    Forward_primer_name = (df['Forward_primer'][sample])[:-len(Forward_primer_number)]
    Reverse_primer_name = (df['Reverse_primer'][sample])[:-len(Reverse_primer_number)]
    
    # Get barcode sequences
    df.at[sample,'ForwardBarcode'] = primer_sequence_library[Forward_primer_name][primer_sequence_library[Forward_primer_name]['Forward_primer']==Forward]['Barcode_Forward_Primer'].to_string(index = False)
    df.at[sample,'ReverseBarcode'] = primer_sequence_library[Reverse_primer_name][primer_sequence_library[Reverse_primer_name]['Reverse_primer']==Reverse]['Barcode_Reverse_Primer'].to_string(index = False)
    
    # check if barcoded primer is barcoded M13
    if "M13" in Forward_primer_name:
        Forward_primer_name = ReadMe.loc[ReadMe['Project_info'] == 'Forward_Primer', 'Fill_In'].iloc[0]     
        Forward_primer_sequence = primer_sequence_library['M13_tailed_primers'].loc[primer_sequence_library['M13_tailed_primers']['Primer'] == Forward_primer_name, 'Primer_sequence'].iloc[0]
        if not M13fwd: 
            fwd_tail = primer_sequence_library['M13_tailed_primers'].loc[primer_sequence_library['M13_tailed_primers']['Primer'] == Forward_primer_name, 'Tail_sequence'].iloc[0]
            M13fwd = True       
    else:
        Forward_primer_sequence = primer_sequence_library[Forward_primer_name]['ForwardPrimer'][0]
    if "M13" in Reverse_primer_name:
        Reverse_primer_name = ReadMe.loc[ReadMe['Project_info'] == 'Reverse_Primer', 'Fill_In'].iloc[0]
        Reverse_primer_sequence = primer_sequence_library['M13_tailed_primers'].loc[primer_sequence_library['M13_tailed_primers']['Primer'] == Reverse_primer_name, 'Primer_sequence'].iloc[0]
        if not M13rev:
            rev_tail = primer_sequence_library['M13_tailed_primers'].loc[primer_sequence_library['M13_tailed_primers']['Primer'] == Reverse_primer_name, 'Tail_sequence'].iloc[0]
            M13rev = True
    else:
        Reverse_primer_sequence = primer_sequence_library[Reverse_primer_name]['ReversePrimer'][0]
    # Add primer name to dataframe
    df.at[sample,'ForwardPrimerName'] = Forward_primer_name
    df.at[sample,'ReversePrimerName'] = Reverse_primer_name
    
    # Add primer and barcode sequences to dataframe
    df.at[sample,'ForwardPrimerSequence'] = Forward_primer_sequence
    df.at[sample,'ReversePrimerSequence'] = Reverse_primer_sequence  

#### Generate SampleIDs
df['#SampleID'] = (
    NIOZnumber + '.' + 
    df ['Forward_primer_number'] + '.' + 
    df ['Reverse_primer_number'])       

#### Add M13 tail sequence to ReadMe
if M13fwd:
    forward_primer_index = ReadMe.loc[(ReadMe == 'Forward_Primer').any(axis=1)].index[0]
    ReadMe = pd.DataFrame(np.insert(ReadMe.values, forward_primer_index+1, values=['M13fwdTail', fwd_tail,'',''], axis=0),columns = ReadMe.columns)

if M13rev:
    reverse_primer_index = ReadMe.loc[(ReadMe == 'Reverse_Primer').any(axis=1)].index[0]
    ReadMe = pd.DataFrame(np.insert(ReadMe.values, reverse_primer_index+1, values=['M13revTail', rev_tail,'',''], axis=0),columns = ReadMe.columns)
    

#### Assemble final mappingfile
# Make an empty dataframe
mf = pd.DataFrame()
# First column: SampleID
mf['#SampleID'] = df['#SampleID']
# Next column: Forward and reverse barcode sequence
mf['ForwardBarcode'] = df['ForwardBarcode']
mf['ReverseBarcode'] = df['ReverseBarcode']
# Next columna: primer names
mf['ForwardPrimerName'] = df['ForwardPrimerName']
mf['ReversePrimerName'] = df['ReversePrimerName']
# Next column: 
mf['ForwardPrimerSequence'] = df['ForwardPrimerSequence']
mf['ReversePrimerSequence'] = df['ReversePrimerSequence']
#### Add metadata from sample_file
for column in sample_file.columns[2:]:
    mf[column] = sample_file[column]

## Save mapping_file
# Save file as RUNID_mapping_file.txt, tab delimited and without the index
mf.to_csv(folder_path + NIOZnumber + "_mapping_file.txt", sep="\t", index=False)

# Saves the ProjectInfo sheet as a .txt file at the same location as the mapping file
ReadMe.to_csv(folder_path + NIOZnumber + "_README!.txt", sep="\t", index=False)

print(f"Your mappingfile can be found here: \n{folder_path}")