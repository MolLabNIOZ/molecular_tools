import pandas as pd


## Import excel sheet
volume_data = pd.read_excel(
    "//zeus.nioz.nl/mmb/molecular_ecology/mollab_team/Projects/2025/MMB/Helge/Gaia/NIOZ418_equimolar_pooling_results_second_try.xlsx",
    sheet_name='part 2')

# Make list of volumes (rounded to 2 decimals)
dna = round(volume_data['water_volume'],2).to_list()

# remove NaNs
dna_noNaNs = [x for x in dna if pd.isnull(x) == False and x != 'nan']

print(dna_noNaNs)

print("\nIs it correct that you have ",len(dna_noNaNs), " samples?")

print("\nCopy_Paste the list with volumes to your OT2 protocol")