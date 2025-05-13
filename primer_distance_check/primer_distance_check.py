# -*- coding: utf-8 -*-
"""
Check a primer set for distance between barcodes
Levenshtein.distance Calculates the minimum number of insertions, deletions, 
and substitutions required to change one sequence into the other according to 
Levenshtein with custom costs for insertion, deletion and substitution.
"""
import Levenshtein as LT
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 

# Make an empty dict and df
barcodes = {}

# Open the file
with open('C:/Users/mbrouwer/OneDrive - NIOZ/Documenten/GitHub/molecular_tools/mapping_file_creator/LCO_HCO_barcodes.txt') as primer_info:
    # Loop through lines of the file
    for i, line in enumerate(primer_info):
        # Skip the header line
        if not i == 0:
            # Get primer number and barcode and append to dict
            fwd_primer_number = line.split()[0][-3:]
            fwd_primer_barcode = line.split()[1]
            barcodes[fwd_primer_number] = fwd_primer_barcode
            
            rev_primer_number = line.split()[3][-3:]
            rev_primer_barcode = line.split()[4]
            barcodes[rev_primer_number] = rev_primer_barcode

# loop through primers and make a list
primers = []
for primer in barcodes:
    primers.append(primer)
# Make an array to compare primers
barcode_array = pd.DataFrame(columns = primers)
barcode_array.insert(0, 'primer', primers)
barcode_array.set_index('primer', inplace=True)


# loop through primers
for primer_1 in barcodes:
    # get barcode sequence
    barcode_1 = barcodes[primer_1]
    
    distances = []
    # loop through primers again
    for primer_2 in barcodes:
        # get barcode sequence to compare with
        barcode_2 = barcodes[primer_2]
        
        # Check minimum number of insertions, deletions, and substitutions required to change one sequence into the other
        distance = LT.distance(barcode_1, barcode_2)
        distances.append(distance)
    
    # Add distances to array
    barcode_array[primer_1] = distances

# # For manual search
# # Sum up distances, to see which have the highest distances
# barcode_array['sum_distances'] = barcode_array.sum(axis = 1)
# # Sort basted on that sum
# barcode_array = barcode_array.sort_values('sum_distances', ascending=False)

    

# # Get an ok set
def ok_primer_set(barcode_array, number_of_primers, max_distance, checkfwdrev):
    ok_primer_set = []
    while len(ok_primer_set) < number_of_primers:
        # Loop through rows of array
        for primer_1, distances in barcode_array.iterrows():
            # Loop through columns of that row
            for primer_2, distance in distances.iteritems():
                if distance >= max_distance:
                    print(primer_2)
    
    
    return ok_primer_set
    
    
ok_primer_set = ok_primer_set(barcode_array = barcode_array, number_of_primers = 10, max_distance = 6, checkfwdrev = False)

        
        
    




    


        
    

        
    
    
        
    
    
        
        
        
    
    


