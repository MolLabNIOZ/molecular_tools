# =============================================================================
# Author:   Maartje Brouwer
# Goal:     Automating mapping file generation
# Date:     210318
# =============================================================================

## Uptil now we had to generate mapping files for every sequencing run by hand.
## Because we always make unique primer-pair combinations for every run, this
## needs to be done for every run seperately.
## In excel this is a lot of work, and often goes wrong. Therefore I would like
## to automate this process.

#### Import needed packages
import pandas as pd       # to be able to work with tables
from Bio.Seq import Seq   # to be able to do compl_rev

# !!! Set variables for your mappingfile
sample_file = 'NIOZ197.csv'
  # Sample_file must be .csv and column 1 must be SampleID
  # SampleID should be NIOZ#.fwprimer#.rvprimer# --> e.g. NIOZ197.001.252
  # After SampleID, you can add as many columns of metadata as you like
  # This metadata will be added to the end of the mappingfile
fw_primer = '515F'
rv_primer = '926RBC'    # '806RB' is also available
  # For any other primer we need to make a new .csv list
  # primer list must look like this: 3 columns:
  # $1 Forward_primer or Reverse_primer --> primername _Golay### 
  # $2 = ForwardPrimer or ReversePrimer --> primer sequence
  # $3 = LinkerPrimerSequence or Barcode_Reverse_Primer --> barcode sequence
  # i.e.: 
  # Reverse_primer;ReversePrimer;Barcode_Reverse_Primer
  # 926RBC_Golay002;CCGYCAATTYMTTTRAGTTT;CTTCCAACTCAT

  
# Import needed files
sample_list = pd.read_csv(sample_file, delimiter=';')
fw_primers = pd.read_csv(fw_primer + '.csv', delimiter=';')
rv_primers = pd.read_csv(rv_primer + '.csv', delimiter=';')

## Get the primer pair numbers from the sampleID
# make an empty dataframe with the desired columns
primerIDs = pd.DataFrame(columns=['SampleID', 'Forward_primer', 'Reverse_primer'])
# append the sampleIDs from the sample_list
primerIDs ['SampleID'] = sample_list['SampleID']
# split fwd and rev primer numbers and add to right column
primerIDs ['Forward_primer'] = (
    (sample_list['SampleID'].str.split('.',expand=True))[1].astype(str))
primerIDs ['Reverse_primer'] = (
    (sample_list['SampleID'].str.split('.',expand=True))[2].astype(str))
# add primer names to the primer numbers
primerIDs ['Forward_primer'] = fw_primer + '_Golay' + primerIDs ['Forward_primer']
primerIDs ['Reverse_primer'] = rv_primer + '_Golay' + primerIDs ['Reverse_primer']

## Add primer sequence and barcode sequences to primerIDs
primerIDs = pd.merge(primerIDs, fw_primers)
primerIDs = pd.merge(primerIDs, rv_primers)

## Get complement reverse of reverse primer barcode
# Make empty new column
primerIDs ['RevComplReverseBarcodesequence'] = ''
# For every sample get the rev_compl of the reverse primer barcode
for sample in primerIDs.index:
    Barcode_Reverse_Primer = Seq(primerIDs['Barcode_Reverse_Primer'][sample])
    compl_rev_barcode = Barcode_Reverse_Primer.reverse_complement()
    # add reverse complement to the dataframe
    primerIDs.at[sample,'RevComplReverseBarcodesequence'] = compl_rev_barcode

## Construct BarcodeSequence (barcode_fwd + revcompl_barcode_rev)
# Make empty new column
primerIDs ['BarcodeSequence'] = ''
# For every sample make the BarcodeSequence
for sample in primerIDs.index:
    BarcodeSequence = ((primerIDs['Barcode_Forward_Primer'][sample]) + 
                        (primerIDs['RevComplReverseBarcodesequence'][sample]))
    # add BarcodeSequence to the dataframe
    primerIDs.at[sample,'BarcodeSequence'] = BarcodeSequence

## Construct InputFileName (full_name_fwd + full_name_rev + .fasta)
# Make empty new column
primerIDs ['InputFileName'] = ''
# For every sample make the InputFileName
for sample in primerIDs.index:
    Forward_Name = (primerIDs['Forward_primer'][sample])
    Reverse_Name = (primerIDs['Reverse_primer'][sample])
    InputFileName = Forward_Name + '.' + Reverse_Name + '.fasta'
    # add BarcodeSequence to the dataframe
    primerIDs.at[sample,'InputFileName'] = InputFileName

## Assemble final mappingfile
# Make an empty dataframe
mf = pd.DataFrame()
# Column 1: SampleID
mf.insert(0, 'SampleID', sample_list['SampleID'])
# Column 2: BarcodeSequence
mf.insert(1, 'BarcodeSequence', primerIDs['BarcodeSequence'])
# Column 3, 4, 5, 6, 7:
# LinkerPrimerSequence, Barcode_Forward_Primer, Barcode_Reverse_Primer, 
# ReversePrimer, RevComplReverseBarcodesequence
mf.insert(2, 'LinkerPrimerSequence', primerIDs['LinkerPrimerSequence'])
mf.insert(3, 'Barcode_Forward_Primer', primerIDs['Barcode_Forward_Primer'])
mf.insert(4, 'Barcode_Reverse_Primer', primerIDs['Barcode_Reverse_Primer'])
mf.insert(5, 'ReversePrimer', primerIDs['ReversePrimer'])
mf.insert(6, 
          'RevComplReverseBarcodesequence', 
          primerIDs['RevComplReverseBarcodesequence'])
# Column 6: InputFileName
mf.insert(7, 'InputFileName', primerIDs['InputFileName'])
# Column 7 and on... extra data from sample_list
for column in sample_list.columns[1:]:
    mf[column] = sample_list[column]

## Save mapping_file
# Extraxt NIOZ number from SampleID
RUNID =(sample_list['SampleID'].str.split('.',expand=True))[0].astype(str).values[0]
# Save file as RUNID_mapping_file.txt, tab delimited and without the index
mf.to_csv(RUNID + "_mapping_file.txt", sep="\t", index=False)