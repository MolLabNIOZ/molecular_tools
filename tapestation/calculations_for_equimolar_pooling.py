"""
This script uses tapestation results for calculating values for equimolar 
pooling
!!! This only works if the configuration of the tapestationplate is equal to the 
PCR plate

Input should be a compactRegionTable as produced by the tapestation (.csv) and
the total amount of DNA (ng) you want in your final pool.

Output can either be the recommendation to dilute the samples first, or to
directly pool the undiluted PCR products.
Furthermore, output will be a list with volumes to pool and if necesarry
the volumes of water/DNA needed to dilute the samples.
"""

# Import needed packages=======================================================
# For working with dataframes we need pandas
import pandas as pd
# To be able to exit the script when the samples are not sufficient we need sys
from sys import exit
# To do some math stuff such ass rounding up etc we need math
import math
# =============================================================================

# Variables to set ============================================================
#### Where is the compactRegionTable .csv located?
filepath = 'C:/Users/rdebeer/Downloads/50_xc-b-jh_240219_NIOZ373 - 2024-02-19 - 12-08-56-D1000-PLATE1-D1000_compactRegionTable.csv'

#### How much PCR product is available (µL)
PCR_volume = 15

#### How much DNA do you want to send for sequencing? (ng)
total_ng = 800
    # The script multiplies this by 2, to take into account you will loose DNA
    # during clean-up

#### If necesarry, how many samples would you dilute by hand, before making an
  ## entire new plate?
max_dilutions_by_hand = 10
# =============================================================================

# Data analysis ===============================================================
#### Read the compactRegionTable .csv and put into a dataframe
data = pd.read_csv(filepath, encoding='unicode-escape')

#### Get a list with all concentrations
concentrations = data['Conc. [ng/µl]'].tolist()
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
        if concentration * PCR_volume < (total_ng * 2) / number_of_samples:
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
ng_per_sample = (total_ng * 2) / len(concentrations)
# Check if a sufficient amount of samples contribute to the pool (>50%)
if len(concentrations) < original_number_of_samples / 2:
    print("Less than halve of your samples has a sufficient amount of DNA to a"
          "dd to the pool. I suggest you either choose a lower [total_ng], add"
          " some samples to your sequencing lane, or re-PCR your samples to ge"
          "t a larger volume.")
    exit()


#### Calculate how many samples need diluting
# you want to pool at least 10 µl of each sample, to make it most accurate
data['dilution_ratio'] = ''
samples_to_dilute = []
preferred_max_concentration = ng_per_sample / 10
# Based on the preferred_max_concentratio check for each sample, whether they 
# need to be diluted and add the dilution ratio to the dataframe.
for sample in data.index:
    # For each sample check the current concentration
    concentration = data['Conc. [ng/µl]'][sample]
    # Calculate the dilution ratio, based on the earlier determined max conc.
    dilution_ratio = concentration / preferred_max_concentration
    # Only if dilution_ratio > 1 a sample needs to be diluted
    if dilution_ratio > 1:
        data.at[sample,'dilution_ratio'] = dilution_ratio
        samples_to_dilute.append(sample)
# How many samples need to be diluted?
number_of_samples_to_dilute = len(samples_to_dilute)
# Recommendation to do it by hand or robot, based on max_dilutions_by_hand
if number_of_samples_to_dilute > max_dilutions_by_hand:
    print (
           str(number_of_samples_to_dilute) + 
           " samples need to be diluted before pooling. "
           "I suggest you let M-O do this for you. This script provides you wi"
           "th a list of volumes for water and a list of volumes for samples. "
           "You can insert these in sample_dilution.py."
           "\nThis will result in an entire new 96-wells plate containing for "
           "each sample either the original sample (these will be transferred "
           "without diluting) or a dilution. This way you can use the entire n"
           "ew 'dilution' plate for equimolar pooling.")
elif number_of_samples_to_dilute == 0:
    print ("No samples need to be diluted. You can use the original PCR produc"
           "ts for equimolar pooling.")
else:
    print (
           str(number_of_samples_to_dilute) + 
           " samples need to be diluted before pooling. "
           "I suggest you do this by hand. " 
           "This script provides you with a list of volumes for water and a li"
           "st of volumes for samples. Make these dulutions in PCR strips. In "
           " the volume list this script provides next, the wells where sample"
           "s that need diluting are located will be skipped and the dilutions"
           " will be added at the end.")

#### Add a column with DNA and water volumes for diluting
# Make 2 new empty columns
data['DNA_volume'] = ''
data['water_volume'] = ''
# If samples need to be diluted, provide the info how to do that
if number_of_samples_to_dilute > 0:
    for sample in data.index:
        dilution_ratio = data['dilution_ratio'][sample]
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
                          f"{data['Sample Description'][sample]}, located in "
                          f"well {data['WellId'][sample]} is too high."
                          f" I suggest you dilute this in the original PCR "
                          f"plate first, adjust the concentration accordingly "
                          f" and run this script again.")
                    exit()
                    
                          
        else:
            DNA_volume = PCR_volume
            water_volume = 0
        # Add values for water and DNA to the dataframe
        data.at[sample,'DNA_volume'] = DNA_volume
        data.at[sample,'water_volume'] = float("%.2f" % water_volume) #2 decimals
        
    #### Print the water and DNA volume lists for diluting. 
    if number_of_samples_to_dilute > max_dilutions_by_hand:
        # M-O will do the dilutions, the lists can be pasted into 
        # sample_dilution.py.
        print("\nNext you will see 2 lists. The first list contains all "
              "volumes of water needed for diluting the samples. The second "
              "list contains volumes of PCR product that need to be "
              "transferred to the new 'dilution' plate. Copy these lists and "
              "paste them into sample_dilution.py"
              "\n\nWater_volumes:")
        print(data['water_volume'].tolist())
        print("\nDNA_volumes:")
        print(data['DNA_volume'].tolist())
    else:
        # These should be diluted by hand, in PCR strips
        print("\n"
              "The following samples need to be diluted. Please do this by "
              "hand.")
        for sample in data.index:
            dilution_ratio = data['dilution_ratio'][sample]
            if (isinstance(dilution_ratio, int) 
                or 
                isinstance(dilution_ratio, float)):
                print("sample " + str(data['Sample Description'][sample]) + 
                      " needs to be diluted. This sample is located in well " +
                      str(data['WellId'][sample]) + " of the PCR_plate. Add to"
                      " a PCR_strip_tube " + str(data['water_volume'][sample]) 
                      + "µL of water and " + str(data['DNA_volume'][sample]) + 
                      "µL of PCR product\n")
                        
# If no samples need to be diluted, DNA_volume = PCR_volume and water_volume = 0
else:
    for sample in data.index:
        data.at[sample,'DNA_volume'] = PCR_volume
        data.at[sample,'water_volume'] = 0

#### Add columns with volumes to use for equimolar pooling and with extra info
# And some additional info to put in the mapping_file
# Make new empty columns
data['final_concentration'] = ''
data['µL_pooled'] = ''
data['ng_pooled'] = ''
data['diluted_before_pooling'] = ''
data[''] = ''
data['pool_information'] = ''
data['values'] = ''
# For every sample calculate the concentration after dilution, 
# calculations work for both diluted and undiluted samples in one batch
if (number_of_samples_to_dilute == 0 
    or 
    number_of_samples_to_dilute > max_dilutions_by_hand):
    for sample in data.index:
        original_concentration = data['Conc. [ng/µl]'][sample]
        DNA_volume = data['DNA_volume'][sample]
        water_volume = data['water_volume'][sample]
        
        # Add final_concentrations to the df
        final_concentration = (original_concentration * DNA_volume) / (water_volume + DNA_volume)
        data.at[sample,'final_concentration'] = final_concentration
        
        # Add volumes to pool to the df
        µL_pooled = float("%.2f" % (ng_per_sample / final_concentration))
        if µL_pooled > PCR_volume:
            µL_pooled = PCR_volume
        data.at[sample,'µL_pooled'] = µL_pooled
        
        # Add total ng per sample pooled to the df
        ng_pooled = final_concentration * µL_pooled
        data.at[sample,'ng_pooled'] = ng_pooled        
        
        # Add info whether a sample has been diluted to the df
        if water_volume > 0:
            diluted = True
        else:
            diluted = False
        data.at[sample,'diluted_before_pooling'] = diluted

# If you did some dilutions by hand, the originals will be skipped and you'll 
# get some separate volumes to pool by hand.
else:
    extra_pool = []
    for sample in data.index:
        dilution_ratio = data['dilution_ratio'][sample]
        original_concentration = data['Conc. [ng/µl]'][sample]
        DNA_volume = data['DNA_volume'][sample]
        water_volume = data['water_volume'][sample]
        final_concentration = (original_concentration * DNA_volume) / (water_volume + DNA_volume)
        µL_pooled = float("%.2f" % (ng_per_sample / final_concentration))
        if µL_pooled > PCR_volume:
            µL_pooled = PCR_volume
            dilution_ratio = data['dilution_ratio'][sample]
        if (isinstance(dilution_ratio, int) 
            or 
            isinstance(dilution_ratio, float)):
            
            µL_pooled = 0
        data.at[sample,'final_concentration'] = final_concentration
        data.at[sample,'µL_pooled'] = µL_pooled

# Define a dictionary to map pool information to their corresponding values
data_mappings = {
    'total µl pooled': data['µL_pooled'].sum(),  # Calculate the total µL pooled
    'total ng pooled': data['ng_pooled'].sum(),  # Calculate the total ng pooled
    'total pb buffer needed (µl)': data['µL_pooled'].sum() * 5,  # Calculate the total pb buffer needed in µl
    'total ph indicator needed (µl)': (data['µL_pooled'].sum() * 5) / 250  # Calculate the total ph indicator needed in µl
}

# Loop over the items in the dictionary and assign values to the dataframe
for i, (info, value) in enumerate(data_mappings.items()):
    data.at[i, 'pool_information'] = info  # Assign pool information to the 'pool_information' column
    data.at[i, 'values'] = value  # Assign values to the 'values' column

#### Print the list with volumes to pool. 
print("\n"
      "Following is a list with volumes you can use for equimolar pooling. Cop"
      "y this list and paste it in equimolar_pooling.py \n\nµL to pool:")
print(data['µL_pooled'].tolist())

#### Save the calculations in a new csv file
filepath = '/'.join(filepath.split('/')[:-1]) + '/'
data.to_csv(filepath + 'equimolar_pooling_results.csv', index=False)


    