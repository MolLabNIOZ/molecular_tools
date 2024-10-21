import pandas as pd


## Import excel sheet
volume_data = pd.read_excel(
    "//zeus.nioz.nl/mmb/molecular_ecology/mollab_team/Projects/2024/COS/JvG/Roos/NIOZ399_equimolar_pooling_results.xlsx",
    sheet_name='NIOZ399_equimolar_pooling_resul')

# Make list of volumes (rounded to 2 decimals)
dna = round(volume_data['50ng'],2).to_list()

# remove NaNs
dna_noNaNs = [x for x in dna if pd.isnull(x) == False and x != 'nan']

print(dna_noNaNs)

print("\nIs it correct that you have ",len(dna_noNaNs), " samples?")

print("\nCopy_Paste the list with volumes to your OT2 protocol")