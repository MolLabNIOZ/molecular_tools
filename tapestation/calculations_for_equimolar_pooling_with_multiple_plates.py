"""
This script uses tapestation results for calculating values for equimolar 
pooling
!!! This only works if the configuration of the tapestationplate is equal to the 
PCR plate

Input should be a compactRegionTable as produced by the tapestation (.csv) and
the total amount of DNA (ng) you want in your final pool.

Output can either be the recommendation to dilute the samples first, or to
directly pool the undiluted PCR products.
Furthermore, output will be opentrons protocols for diluting and/or equimolar
pooling
"""

# Variables to set ============================================================
#### What is the folder where your compactRegionTable.csv's are located?
folderpath = 'molecular_tools/tests_RdB/NIOZ999_equimolar_pooling test'

#### How much PCR product is available (µL)
PCR_volume = 25

#### How much DNA do you want to send for sequencing? (ng)
total_ng = 500

#### If necesarry, how many samples would you dilute by hand, before making an
  ## entire new plate?
max_dilutions_by_hand = 10
# =============================================================================

# Import needed packages=======================================================
# For working with dataframes we need pandas
import pandas as pd
# To be able to exit the script when the samples are not sufficient we need sys
from sys import exit
# To do some math stuff such ass rounding up etc we need math
import math
# Imports shutil, this is used for creating a copy of the template file
import shutil
import glob
import os
# =============================================================================

# Data analysis ===============================================================
# Get the NIOZ number from the filename
if 'NIOZ' in folderpath:
    NIOZ_number = 'NIOZ' + folderpath.split('NIOZ',1)[1][:3]
else:
    raise Exception("Make sure your NIOZnumber is in the data filename")

# Loops through the folder, uses glob to get a list of all CSV file paths in the specified folder 
file_names = glob.glob(os.path.join(folderpath,'*.csv'))
# Filter the list of file paths to include only those whose file names contain 'compactRegionTable'
# The `os.path.basename(file)` extracts the file name from the full file path
file_names = [file for file in file_names if 'compactRegionTable' in os.path.basename(file)]


#### Read the compactRegionTable .csv and put into a dict
data_dict = {}
for file in file_names:
    plate_number = 'plate_' + file.split('plate_')[1][:1]
    quantification_file = pd.read_csv(file, encoding = 'unicode-escape')
    quantification_file['plate_number'] = plate_number
    data_dict[plate_number] = quantification_file

# Creates one dataframe from the dictonary created above
data_df = pd.concat(data_dict.values(), ignore_index = True)

#### Defines new map on ZEUS
new_folder = '//lab-mmb.nioz.nl/logs/MolLab_robots/Protocol_database/MO/Generated_protocols/' + NIOZ_number

#### Get a list with all concentrations
concentrations = data_df['Conc. [ng/µl]'].tolist()

#### check how many samples there are in the list
original_number_of_samples = len(concentrations)

#### Deduct from the dataset, how much DNA (ng) to pool per sample
# Check based on the sampleset what samples have sufficient DNA to contribute
def remove_insufficient_samples():
    # Make a recursive loop, that keeps looping until only sufficient samples
    # are left in the list
    recursive_loop = False
    for concentration in concentrations:
        number_of_samples = len(concentrations)
        # Calculate for each sample whether there is a sufficient amount of DNA
        if concentration * PCR_volume < (total_ng) / number_of_samples:
            # If not, remove it from the list
            concentrations.remove(concentration)
        # Because now there is a lower amount of samples the amount of DNA
        # needed per sample is higher and the loop has to be repeated    
            recursive_loop = True
        if recursive_loop:
            remove_insufficient_samples()
        # If no samples are removed, the recursive_loop is not turned to True
        # and the loop stops
        else:
            pass
# Calling the previously written function
remove_insufficient_samples()

# Calculate the ng per sample needed
ng_per_sample = (total_ng) / len(concentrations)
# Check if a sufficient amount of samples contribute to the pool (>50%)
if len(concentrations) < original_number_of_samples / 2:
    print("Less than halve of your samples has a sufficient amount of DNA to a"
          "dd to the pool. I suggest you either choose a lower [total_ng], add"
          " some samples to your sequencing lane, or re-PCR your samples to ge"
          "t a larger volume.")
    exit()

#### Calculate how many samples need diluting
# you want to pool at least 10 µl of each sample, to make it most accurate
data_df['dilution_ratio'] = ''
samples_to_dilute = []
preferred_max_concentration = ng_per_sample / 10
# Based on the preferred_max_concentratio check for each sample, whether they 
# need to be diluted and add the dilution ratio to the dataframe.
for sample in data_df.index:
    # For each sample check the current concentration
    concentration = data_df['Conc. [ng/µl]'][sample]
    # Calculate the dilution ratio, based on the earlier determined max conc.
    dilution_ratio = concentration / preferred_max_concentration
    # Only if dilution_ratio > 1 a sample needs to be diluted
    if dilution_ratio > 1:
        data_df.at[sample,'dilution_ratio'] = dilution_ratio
        samples_to_dilute.append(sample)

# How many samples need to be diluted?
number_of_samples_to_dilute = len(samples_to_dilute)
# Recommendation to do it by hand or robot, based on max_dilutions_by_hand
if number_of_samples_to_dilute > max_dilutions_by_hand:
    print (str(number_of_samples_to_dilute) + 
           " samples need to be diluted before pooling. "
           "I suggest you let me do this for you. I have been so kind to "
           "create the two procols you need for this! You can find them in "
           "the M-O folder and contain your NIOZ number. Good luck!"
           "\nUntil next time! \nM-O")
           
elif number_of_samples_to_dilute == 0:
    print ("No samples need to be diluted. You can use the original PCR produc"
           "ts for equimolar pooling. For pooling I have already created a "
           "protocol.You can find them in "
           "the M-O folder and contain your NIOZ number. Good luck!"
           "\nUntil next time! \nM-O ")
else:
    print (
           str(number_of_samples_to_dilute) + 
           " samples need to be diluted before pooling. "
           "I suggest you do this by hand. " 
           "This script provides you with a list of volumes for water and a li"
           "st of volumes for samples. Make these dilutions in PCR strips. In "
           " the volume list this script provides next, the wells where sample"
           "s that need diluting are located will be skipped and the dilutions"
           " will be added at the end.")

#### Add a column with DNA and water volumes for diluting
# Make 2 new empty columns
data_df['DNA_volume'] = ''
data_df['water_volume'] = ''
# If samples need to be diluted, provide the info how to do that
if number_of_samples_to_dilute > 0:
    for sample in data_df.index:
        dilution_ratio = data_df['dilution_ratio'][sample]
        # For every sample determine how much PCR product to use and if it
        # needs to be diluted, how much water is needed for that
        if dilution_ratio:
            dilution_ratio = float(dilution_ratio)
            # In steps of 10µL, check if it yields sufficient amounts of DNA
            for i in range(1,math.ceil(PCR_volume/10)):
                DNA_volume = 10 * i
                if dilution_ratio * DNA_volume > PCR_volume:
                    break
                else:
                    continue
                
            water_volume = DNA_volume * dilution_ratio - DNA_volume
            while water_volume + DNA_volume > 200:
                water_volume = water_volume / 2
                DNA_volume = DNA_volume / 2
                if DNA_volume < 2.5:
                    print(f"\n!!!\nConcentration of sample " 
                          f"{data_df['Sample Description'][sample]}, located in "
                          f"well {data_df['WellId'][sample]} is too high."
                          f" I suggest you dilute this in the original PCR "
                          f"plate first, adjust the concentration accordingly "
                          f" and run this script again.")
                    exit()
                    
                          
        else:
            DNA_volume = PCR_volume
            water_volume = 0
        # Add values for water and DNA to the dataframe
        data_df.at[sample,'DNA_volume'] = DNA_volume
        data_df.at[sample,'water_volume'] = float("%.2f" % water_volume) #2 decimals
# If no samples need to be diluted, DNA_volume = PCR_volume and water_volume = 0
else:
    for sample in data_df.index:
        data_df.at[sample,'DNA_volume'] = PCR_volume
        data_df.at[sample,'water_volume'] = 0

#### Add columns with volumes to use for equimolar pooling and with extra info
# And some additional info to put in the mapping_file
# Make new empty columns
data_df['final_concentration'] = ''
data_df['µL_pooled'] = ''
data_df['ng_pooled'] = ''
data_df['diluted_before_pooling'] = ''
data_df['ng_equimolar'] = ''
data_df[''] = ''
data_df['pool_information'] = ''
data_df['values'] = ''
# For every sample calculate the concentration after dilution, 
# calculations work for both diluted and undiluted samples in one batch
if (number_of_samples_to_dilute == 0 
    or 
    number_of_samples_to_dilute > max_dilutions_by_hand):
    for sample in data_df.index:
        original_concentration = data_df['Conc. [ng/µl]'][sample]
        DNA_volume = data_df['DNA_volume'][sample]
        water_volume = data_df['water_volume'][sample]
        
        # Add final_concentrations to the df
        final_concentration = (original_concentration * DNA_volume) / (water_volume + DNA_volume)
        data_df.at[sample,'final_concentration'] = final_concentration
        
        # Add volumes to pool to the df
        µL_pooled = float("%.2f" % (ng_per_sample / final_concentration))
        if µL_pooled > PCR_volume:
            µL_pooled = PCR_volume
        data_df.at[sample,'µL_pooled'] = µL_pooled
        
        # Add total ng per sample pooled to the df
        ng_pooled = final_concentration * µL_pooled
        data_df.at[sample,'ng_pooled'] = ng_pooled        
        
        # Add info whether a sample has been diluted to the df
        if water_volume > 0:
            diluted = True
        else:
            diluted = False
        data_df.at[sample,'diluted_before_pooling'] = diluted
        
        # Add info whether a sample has sufficient DNA to reach the asked ng
        if math.ceil(ng_pooled) >= ng_per_sample:
            equimolar = True
        else:
            equimolar = False
        data_df.at[sample,'ng_equimolar'] = equimolar
        
# If you did some dilutions by hand, the originals will be skipped and you'll 
# get some separate volumes to pool by hand.
else:
    extra_pool = []
    for sample in data_df.index:
        dilution_ratio = data_df['dilution_ratio'][sample]
        original_concentration = data_df['Conc. [ng/µl]'][sample]
        DNA_volume = data_df['DNA_volume'][sample]
        water_volume = data_df['water_volume'][sample]
        final_concentration = (original_concentration * DNA_volume) / (water_volume + DNA_volume)
        µL_pooled = float("%.2f" % (ng_per_sample / final_concentration))
        if µL_pooled > PCR_volume:
            µL_pooled = PCR_volume
            dilution_ratio = data_df['dilution_ratio'][sample]
        if (isinstance(dilution_ratio, int) 
            or 
            isinstance(dilution_ratio, float)):
            
            µL_pooled = 0
        data_df.at[sample,'final_concentration'] = final_concentration
        data_df.at[sample,'µL_pooled'] = µL_pooled

# Define a dictionary to map pool information to their corresponding values
data_mappings = {
    'total µl pooled': data_df['µL_pooled'].sum(),  # Calculate the total µL pooled
    'total ng pooled': data_df['ng_pooled'].sum(),  # Calculate the total ng pooled
    'total pb buffer needed (µl)': data_df['µL_pooled'].sum() * 5,  # Calculate the total pb buffer needed in µl
    'total ph indicator needed (µl)': (data_df['µL_pooled'].sum() * 5) / 250  # Calculate the total ph indicator needed in µl
    }

# Loop over the items in the dictionary and assign values to the dataframe
for i, (info, value) in enumerate(data_mappings.items()):
    data_df.at[i, 'pool_information'] = info  # Assign pool information to the 'pool_information' column
    data_df.at[i, 'values'] = value  # Assign values to the 'values' column
    
#### Save the calculations in a new csv file
data_df.to_csv(f'{folderpath}/{NIOZ_number}_equimolar_pooling_results.csv', index = False, encoding = 'ansi')

# =============================================================================
### CREATING SEPARETE DATAFRAMES ==============================================
# =============================================================================
# Create unique list of the plates
plates = data_df.plate_number.unique()

# Creates a dataframe dictionary to store the dataframes per plate in
final_data_dict = {elem : pd.DataFrame() for elem in plates}
for key in final_data_dict.keys():
    final_data_dict[key] = data_df[:][data_df.plate_number == key]

# =============================================================================
### CREATING NEW DILUTING PROTOCOLS ===========================================
# =============================================================================
#### Print the water and DNA volume lists for diluting. 
if number_of_samples_to_dilute > max_dilutions_by_hand:
    for dataframe_name, dataframe in final_data_dict.items():
        sample_volumes = (dataframe['DNA_volume'].tolist())
        water_volumes = (dataframe['water_volume'].tolist())
        # Locating the template file
        # The directory for the new file with the name it should get
        template_file = 'OT2/Protocol_database/MO/pool_template_protocols/Sample_dilution_protocol_template.py'
        if not os.path.exists(new_folder):
            os.mkdir(new_folder)
        destination_pathway =  new_folder + '/' + NIOZ_number +  f'_{dataframe_name}_sample_dilution.py'
        # Creates the copy of the right templates
        shutil.copy(template_file, destination_pathway)
        # Replace placeholders with lists of volumes
        search_sample = '<Sample_volumes>'
        search_water = '<Water_volumes>'
        search_NIOZ_number = '<NIOZ_NUMBER>'
        with open (destination_pathway, 'r') as file:
            dilution = file.read()
            replace_water_sample = dilution.replace(search_sample, str(sample_volumes)).replace(search_water, str(water_volumes)).replace(search_NIOZ_number, NIOZ_number) 
        # Write modified content back to the file            
        with open(destination_pathway, 'w') as file:
            file.write(replace_water_sample)
else:
    # These should be diluted by hand, in PCR strips
    print("\n"
          "The following samples need to be diluted. Please do this by "
          "hand.")
    for sample in data_df.index:
        dilution_ratio = data_df['dilution_ratio'][sample]
        if (isinstance(dilution_ratio, int) 
            or 
            isinstance(dilution_ratio, float)):
            print("sample " + str(data_df['Sample Description'][sample]) + 
                  " needs to be diluted. This sample is located in well " +
                  str(data_df['WellId'][sample]) + " of the PCR_plate. Add to"
                  " a PCR_strip_tube " + str(data_df['water_volume'][sample]) 
                  + "µL of water and " + str(data_df['DNA_volume'][sample]) + 
                  "µL of PCR product\n")
# =============================================================================
### CREATING NEW POOLING PROTOCOLS ============================================
# =============================================================================  
DNA_volumes = []
for dataframe_name, dataframe in final_data_dict.items():
    DNA_volume = (dataframe['µL_pooled'].tolist())
    DNA_volumes.append(DNA_volume)
    
#### Make opentrons protocols for pooling
# Locating the template file
# The directory for the new file with the name it should get
pooling_template_file = 'OT2/Protocol_database/MO/pool_template_protocols/Equimolar_pooling_with_multiple_plates.py'
destination_pathway = new_folder + '/' + NIOZ_number + '_equimolar_pooling.py'
# Creates the copy of the right templates
shutil.copy(pooling_template_file, destination_pathway)

search_DNA = '<DNA_volumes>'
search_NIOZ_number = '<NIOZ_NUMBER>'

# Replace placeholders                
with open (destination_pathway, 'r') as file:
    pooling = file.read()
    replace_DNA = pooling.replace(search_DNA, str(DNA_volumes)).replace(search_NIOZ_number, NIOZ_number) 
# Write modified content back to the file            
with open(destination_pathway, 'w') as file:
    file.write(replace_DNA)